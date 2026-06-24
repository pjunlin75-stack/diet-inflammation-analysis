#!/usr/bin/env python3
"""
NHANES Chronic Psychological Stress & Allostatic Load × Mortality Analysis
==========================================================================
Examines whether chronic psychological stress phenotype (CPS, top-decile
PHQ-9) accelerates all-cause mortality via allostatic load pathways —
cardiovascular dysregulation (hypertension) and systemic inflammation
(elevated hs-CRP).

Theoretical framework:
  McEwen (1998) Allostatic Load model — cumulative physiological wear from
  chronic stress exposure, measured through cardiovascular and inflammatory
  biomarkers as secondary mediators of the stress–mortality pathway.

Data sources (NHANES, publicly available):
  - DPQ  : Patient Health Questionnaire (PHQ-9) depression screener
  - BPX  : Blood pressure examination
  - CRP  : C-reactive protein (standard or high-sensitivity)
  - DEMO : Demographics (age, sex)
  - Linked Mortality : NCHS public-use linked mortality through 2019

NOTE: 2011-2014 cycles excluded by default — CDC did not release standalone
      CRP lab files for those years.

Author : Luster Dawn
License: MIT
"""

from __future__ import annotations

import io
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from lifelines import CoxPHFitter, KaplanMeierFitter

warnings.filterwarnings("ignore", category=FutureWarning)


# ─────────────────────────────────────────────────────────────────────
# §1  Configuration
# ─────────────────────────────────────────────────────────────────────

@dataclass
class NHANESConfig:
    """Centralised knobs — change cycles/thresholds here, nowhere else."""

    # Cycles with BOTH PHQ-9 AND standalone CRP data
    # 2011-2014 excluded: no CRP lab file those years
    cycles: list[str] = field(default_factory=lambda: [
        "2005-2006", "2007-2008", "2009-2010",
        "2015-2016", "2017-2018",
    ])

    # PHQ-9 item variable names (DPQ010–DPQ090)
    phq9_items: list[str] = field(default_factory=lambda: [
        f"DPQ0{i}0" for i in range(1, 10)
    ])

    # CPS quantile cutoff (top 10 % = chronic psychological stress phenotype)
    quantile_cutoff: float = 0.90

    # CRP threshold for elevated systemic inflammation (mg/L)
    # Per AHA/CDC guideline: ≥3.0 mg/L = high cardiovascular risk
    crp_high_threshold: float = 3.0

    # Hypertension thresholds (mmHg)
    sbp_hypertension: float = 140.0
    dbp_hypertension: float = 90.0

    # Output directory
    output_dir: Path = Path("output")

    # CDC URLs
    xpt_base: str = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"
    mortality_base: str = (
        "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/"
        "datalinkage/linked_mortality"
    )

    def __post_init__(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────
# §2  Data Loader
# ─────────────────────────────────────────────────────────────────────

class NHANESLoader:
    """Downloads NHANES SAS transport (.XPT) and mortality (.dat) files."""

    # cycle → (suffix letter, start year, CRP filename prefix)
    CYCLE_META: dict[str, dict] = {
        "2005-2006": {"sfx": "D", "year": "2005", "crp_prefix": "CRP"},
        "2007-2008": {"sfx": "E", "year": "2007", "crp_prefix": "CRP"},
        "2009-2010": {"sfx": "F", "year": "2009", "crp_prefix": "CRP"},
        # 2011-2014: no standalone CRP file
        "2015-2016": {"sfx": "I", "year": "2015", "crp_prefix": "HSCRP"},
        "2017-2018": {"sfx": "J", "year": "2017", "crp_prefix": "HSCRP"},
    }

    REQUEST_TIMEOUT = 120  # seconds

    def __init__(self, cfg: NHANESConfig):
        self.cfg = cfg

    # ── public ──────────────────────────────────────────────────────

    def load_all(self) -> pd.DataFrame:
        """Pool requested cycles into one DataFrame per component,
        inner-join on SEQN, return the master table."""
        frames = {k: [] for k in ("dpq", "bpx", "crp", "demo")}

        for cycle in self.cfg.cycles:
            meta = self.CYCLE_META[cycle]
            sfx, yr = meta["sfx"], meta["year"]
            crp_pfx = meta["crp_prefix"]

            frames["dpq"].append(self._fetch_xpt(yr, f"DPQ_{sfx}"))
            frames["bpx"].append(self._fetch_xpt(yr, f"BPX_{sfx}"))
            frames["crp"].append(self._fetch_xpt(yr, f"{crp_pfx}_{sfx}"))
            frames["demo"].append(self._fetch_xpt(yr, f"DEMO_{sfx}"))

        dpq = pd.concat(frames["dpq"], ignore_index=True)
        bpx = pd.concat(frames["bpx"], ignore_index=True)
        crp = pd.concat(frames["crp"], ignore_index=True)
        demo = pd.concat(frames["demo"], ignore_index=True)
        mort = self._load_mortality()

        # Inner join on SEQN
        master = (
            dpq.merge(bpx, on="SEQN", how="inner")
               .merge(crp, on="SEQN", how="inner")
               .merge(demo, on="SEQN", how="inner")
               .merge(mort, on="SEQN", how="inner")
        )
        print(f"[Loader] Pooled master table: {master.shape[0]:,} rows × "
              f"{master.shape[1]} cols")
        return master

    # ── private ─────────────────────────────────────────────────────

    def _fetch_xpt(self, start_year: str, filename: str) -> pd.DataFrame:
        """Download XPT via the direct-file URL pattern."""
        url = (f"{self.cfg.xpt_base}/{start_year}"
               f"/DataFiles/{filename}.XPT")
        print(f"  ↓ {url}")
        resp = requests.get(url, timeout=self.REQUEST_TIMEOUT)
        resp.raise_for_status()
        buf = io.BytesIO(resp.content)
        df = pd.read_sas(buf, format="xport", encoding="utf-8")
        # Fix SAS XPT encoding: 0 values are sometimes stored as tiny
        # positive floats (~5.4e-79). Round them back to 0.
        numeric_cols = df.select_dtypes(include="number").columns
        for col in numeric_cols:
            mask = (df[col].notna()) & (df[col] != 0) & (df[col].abs() < 1e-70)
            df.loc[mask, col] = 0.0
        return df

    def _load_mortality(self) -> pd.DataFrame:
        """Load per-cycle NCHS linked mortality files.

        Column positions from the official NCHS R program:
        R_ReadInProgramAllSurveys.R on the CDC FTP.
        """
        colspecs = [
            (0, 6),      # SEQN
            (14, 15),    # ELIGSTAT
            (15, 16),    # MORTSTAT
            (16, 19),    # UCOD_LEADING
            (19, 20),    # DIABETES
            (20, 21),    # HYPERTEN
            (42, 45),    # PERMTH_INT
            (45, 48),    # PERMTH_EXM
        ]
        names = [
            "SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING",
            "DIABETES_MORT", "HYPERTEN_MORT",
            "PERMTH_INT", "PERMTH_EXM",
        ]

        # Map cycle → mortality filename
        cycle_to_file = {
            "2005-2006": "NHANES_2005_2006_MORT_2019_PUBLIC.dat",
            "2007-2008": "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
            "2009-2010": "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
            "2015-2016": "NHANES_2015_2016_MORT_2019_PUBLIC.dat",
            "2017-2018": "NHANES_2017_2018_MORT_2019_PUBLIC.dat",
        }

        mort_frames = []
        for cycle in self.cfg.cycles:
            fname = cycle_to_file[cycle]
            url = f"{self.cfg.mortality_base}/{fname}"
            print(f"  ↓ {url}")

            resp = requests.get(url, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
            buf = io.StringIO(resp.text)

            df = pd.read_fwf(
                buf,
                colspecs=colspecs,
                names=names,
                na_values=[".", ""],
            )
            mort_frames.append(df)

        mort = pd.concat(mort_frames, ignore_index=True)

        # Keep only mortality-eligible respondents
        mort = mort[mort["ELIGSTAT"] == 1].copy()
        mort["PERMTH_INT"] = pd.to_numeric(mort["PERMTH_INT"], errors="coerce")
        mort["SURV_YEARS"] = mort["PERMTH_INT"] / 12.0
        mort["MORTSTAT"] = pd.to_numeric(mort["MORTSTAT"], errors="coerce")

        print(f"  [Mortality] {len(mort):,} eligible, "
              f"{int(mort['MORTSTAT'].sum()):,} deceased")
        return mort[["SEQN", "MORTSTAT", "SURV_YEARS"]].dropna()


# ─────────────────────────────────────────────────────────────────────
# §3  Data Cleaner
# ─────────────────────────────────────────────────────────────────────

class DataCleaner:
    """Applies domain-aware cleaning rules to the master table."""

    def __init__(self, cfg: NHANESConfig):
        self.cfg = cfg

    def clean(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()

        # 3a. PHQ-9: recode 7/9 (refused/don't know) → NaN, then sum
        phq_cols = [c for c in self.cfg.phq9_items if c in df.columns]
        for item in phq_cols:
            df[item] = pd.to_numeric(df[item], errors="coerce")
            df.loc[df[item].isin([7, 9]), item] = np.nan

        df["PHQ9_TOTAL"] = df[phq_cols].sum(axis=1, min_count=len(phq_cols))

        # 3b. Blood pressure: average of available readings
        # Standard naming (2005-2016): BPXSY1-4, BPXDI1-4
        # 2017-2018 naming: BPXOSY1-3, BPXODI1-3
        sbp_cols = sorted([c for c in df.columns
                           if c.startswith(("BPXSY", "BPXOSY"))
                           and c[-1].isdigit()])
        dbp_cols = sorted([c for c in df.columns
                           if c.startswith(("BPXDI", "BPXODI"))
                           and c[-1].isdigit()])

        if sbp_cols:
            df["SBP_MEAN"] = (df[sbp_cols]
                              .apply(pd.to_numeric, errors="coerce")
                              .mean(axis=1))
        else:
            df["SBP_MEAN"] = np.nan

        if dbp_cols:
            df["DBP_MEAN"] = (df[dbp_cols]
                              .apply(pd.to_numeric, errors="coerce")
                              .mean(axis=1))
        else:
            df["DBP_MEAN"] = np.nan

        # 3c. CRP — harmonise across naming conventions
        crp_col = None
        for candidate in ("LBXHSCRP", "LBXCRP"):
            if candidate in df.columns:
                crp_col = candidate
                break
        if crp_col:
            df["CRP"] = pd.to_numeric(df[crp_col], errors="coerce")
        else:
            df["CRP"] = np.nan

        # 3d. Demographics
        df["AGE"] = pd.to_numeric(df.get("RIDAGEYR"), errors="coerce")
        df["FEMALE"] = (pd.to_numeric(df.get("RIAGENDR"),
                                      errors="coerce") == 2).astype(float)

        # 3e. Drop rows missing critical fields
        required = ["PHQ9_TOTAL", "SURV_YEARS", "MORTSTAT", "AGE"]
        before = len(df)
        df.dropna(subset=required, inplace=True)
        df = df[df["SURV_YEARS"] > 0].copy()
        print(f"[Cleaner] {before - len(df):,} rows dropped → "
              f"{len(df):,} remain")

        return df


# ─────────────────────────────────────────────────────────────────────
# §4  Feature Engineer (Academic Refinement)
# ─────────────────────────────────────────────────────────────────────

class FeatureEngineer:
    """Derives clinical and epidemiological variables from cleaned data.

    Operationalises the allostatic load framework (McEwen, 1998) using
    NHANES biomarkers as secondary mediators:
      - CPS Phenotype  : chronic psychological stress (PHQ-9 top decile)
      - Hypertension   : cardiovascular allostatic load component
      - Elevated CRP   : inflammatory allostatic load component
      - CPS × Inflam   : synergistic interaction term
    """

    def __init__(self, cfg: NHANESConfig):
        self.cfg = cfg

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()

        # 4a. Chronic Psychological Stress (CPS) Phenotype
        #     — Top decile of PHQ-9 distribution (proxy for high allostatic load)
        cutoff = df["PHQ9_TOTAL"].quantile(self.cfg.quantile_cutoff)
        cutoff = max(cutoff, 1.0)  # Ensure meaningful cutoff
        df["CPS_PHENOTYPE"] = (df["PHQ9_TOTAL"] >= cutoff).astype(int)
        n_cps = df["CPS_PHENOTYPE"].sum()
        print(f"[Features] CPS Phenotype (PHQ-9 ≥ P"
              f"{int(self.cfg.quantile_cutoff * 100)}) "
              f"→ n={n_cps:,} ({100 * n_cps / len(df):.1f}%)")

        # 4b. Allostatic Load Component: Hypertension
        df["HYPERTENSION"] = (
            (df["SBP_MEAN"] >= self.cfg.sbp_hypertension) |
            (df["DBP_MEAN"] >= self.cfg.dbp_hypertension)
        ).astype(int)

        # 4c. Inflammatory Biomarker: Elevated hs-CRP (≥ 3.0 mg/L)
        df["ELEVATED_CRP"] = (
            df["CRP"] >= self.cfg.crp_high_threshold
        ).astype(int)

        # 4d. Synergistic Interaction: CPS × Inflammation
        #     Tests whether chronic stress and systemic inflammation
        #     produce supra-additive mortality risk (allostatic overload)
        df["CPS_x_INFLAM"] = df["CPS_PHENOTYPE"] * df["ELEVATED_CRP"]

        # 4e. Survival outcomes
        df["EVENT"] = df["MORTSTAT"].astype(int)
        df["DURATION"] = df["SURV_YEARS"]

        return df


# ─────────────────────────────────────────────────────────────────────
# §5  Survival Modeler
# ─────────────────────────────────────────────────────────────────────

class SurvivalModeler:
    """Fits Cox proportional-hazards model for allostatic load pathway.

    Covariates:
      CPS_PHENOTYPE  — chronic psychological stress index (PHQ-9 ≥ P90)
      HYPERTENSION   — cardiovascular allostatic load marker
      ELEVATED_CRP   — inflammatory allostatic load marker
      CPS_x_INFLAM   — stress-inflammation synergistic interaction
      AGE, FEMALE     — demographic confounders
    """

    COVARIATES = [
        "CPS_PHENOTYPE",
        "HYPERTENSION",
        "ELEVATED_CRP",
        "CPS_x_INFLAM",
        "AGE",
        "FEMALE",
    ]

    def __init__(self, cfg: NHANESConfig):
        self.cfg = cfg
        self.cph = CoxPHFitter()
        self._fitted = False

    def fit(self, df: pd.DataFrame) -> "SurvivalModeler":
        model_df = df[self.COVARIATES + ["DURATION", "EVENT"]].dropna()
        print(f"[Model] Fitting Cox PH on {len(model_df):,} complete cases")

        self.cph.fit(
            model_df,
            duration_col="DURATION",
            event_col="EVENT",
        )
        self._fitted = True
        return self

    def summary(self) -> pd.DataFrame:
        """Return hazard-ratio table with confidence intervals."""
        assert self._fitted, "Call .fit() first"
        s = self.cph.summary
        s["HR"] = np.exp(s["coef"])
        s["HR_lower"] = np.exp(s["coef lower 95%"])
        s["HR_upper"] = np.exp(s["coef upper 95%"])
        return s[["HR", "HR_lower", "HR_upper", "p"]].round(4)

    def print_summary(self):
        print("\n" + "=" * 65)
        print("  Cox Proportional-Hazards — Hazard Ratio Summary")
        print("=" * 65)
        print(self.summary().to_string())
        print("=" * 65 + "\n")
        self.cph.print_summary(columns=["coef", "exp(coef)", "p"])


# ─────────────────────────────────────────────────────────────────────
# §6  Visualiser
# ─────────────────────────────────────────────────────────────────────

class SurvivalVisualiser:
    """Generates publication-quality survival curves for CPS analysis."""

    PALETTE = {
        "cps_phenotype": "#D7263D",   # warm red
        "control":       "#1B998B",   # teal
    }

    def __init__(self, cfg: NHANESConfig):
        self.cfg = cfg

    def plot_km_curves(self, df: pd.DataFrame, save: bool = True):
        """Kaplan-Meier curves stratified by CPS phenotype."""
        fig, ax = plt.subplots(figsize=(9, 6))
        kmf = KaplanMeierFitter()

        for label, mask, color in [
            ("Control (PHQ-9 < P90)",
             df["CPS_PHENOTYPE"] == 0, self.PALETTE["control"]),
            ("CPS Phenotype (PHQ-9 ≥ P90)",
             df["CPS_PHENOTYPE"] == 1, self.PALETTE["cps_phenotype"]),
        ]:
            subset = df.loc[mask]
            kmf.fit(
                durations=subset["DURATION"],
                event_observed=subset["EVENT"],
                label=label,
            )
            kmf.plot_survival_function(ax=ax, color=color, linewidth=2)

        ax.set_title(
            "Baseline Survival by Chronic Psychological Stress Phenotype\n"
            "(NHANES pooled, NCHS Linked Mortality 2019)",
            fontsize=13, fontweight="bold",
        )
        ax.set_xlabel("Follow-up (years)", fontsize=11)
        ax.set_ylabel("Survival Probability", fontsize=11)
        ax.set_ylim(0, 1.02)
        ax.legend(loc="lower left", fontsize=10, frameon=True)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()

        if save:
            out = self.cfg.output_dir / "km_cps_mortality.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            print(f"[Viz] Saved → {out}")
        plt.show()

    def plot_cox_baseline(self, modeler: SurvivalModeler,
                          save: bool = True):
        """Plot partial effect of CPS phenotype from Cox model."""
        assert modeler._fitted, "Fit model first"

        fig, ax = plt.subplots(figsize=(9, 6))
        modeler.cph.plot_partial_effects_on_outcome(
            covariates="CPS_PHENOTYPE",
            values=[0, 1],
            ax=ax,
            cmap="RdYlGn_r",
        )
        ax.set_title(
            "Cox PH — Partial Effect of CPS Phenotype on Survival",
            fontsize=13, fontweight="bold",
        )
        ax.set_xlabel("Follow-up (years)", fontsize=11)
        ax.set_ylabel("Baseline Survival Probability", fontsize=11)
        ax.legend(["Control", "CPS Phenotype"],
                  loc="lower left", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()

        if save:
            out = self.cfg.output_dir / "cox_partial_cps.png"
            fig.savefig(out, dpi=300, bbox_inches="tight")
            print(f"[Viz] Saved → {out}")
        plt.show()


# ─────────────────────────────────────────────────────────────────────
# §7  Pipeline Orchestrator
# ─────────────────────────────────────────────────────────────────────

class AllostaticLoadPipeline:
    """End-to-end orchestrator — load → clean → feature → model → viz.

    Implements the CPS–allostatic load–mortality pathway analysis.
    """

    def __init__(self, cfg: Optional[NHANESConfig] = None):
        self.cfg = cfg or NHANESConfig()
        self.loader = NHANESLoader(self.cfg)
        self.cleaner = DataCleaner(self.cfg)
        self.engineer = FeatureEngineer(self.cfg)
        self.modeler = SurvivalModeler(self.cfg)
        self.viz = SurvivalVisualiser(self.cfg)

        self.raw: Optional[pd.DataFrame] = None
        self.clean: Optional[pd.DataFrame] = None
        self.analytic: Optional[pd.DataFrame] = None

    def run(self):
        print("\n╔══════════════════════════════════════════════════════════╗")
        print("║  NHANES CPS × Allostatic Load × Mortality Pipeline      ║")
        print("╚══════════════════════════════════════════════════════════╝\n")

        # Step 1 — Load & merge
        print("─── Stage 1: Loading NHANES components ───")
        self.raw = self.loader.load_all()

        # Step 2 — Clean
        print("\n─── Stage 2: Cleaning ───")
        self.clean = self.cleaner.clean(self.raw)

        # Step 3 — Feature engineering
        print("\n─── Stage 3: Feature Engineering ───")
        self.analytic = self.engineer.transform(self.clean)

        # Step 4 — Cox PH model
        print("\n─── Stage 4: Cox PH Modeling ───")
        self.modeler.fit(self.analytic)
        self.modeler.print_summary()

        # Step 5 — Visualisation
        print("\n─── Stage 5: Visualisation ───")
        self.viz.plot_km_curves(self.analytic)
        self.viz.plot_cox_baseline(self.modeler)

        # Save analytic dataset
        out_csv = self.cfg.output_dir / "analytic_dataset.csv"
        self.analytic.to_csv(out_csv, index=False)
        print(f"\n[Pipeline] Analytic dataset saved → {out_csv}")
        print("[Pipeline] ✓ Complete.\n")
        return self


# ─────────────────────────────────────────────────────────────────────
# §8  Entry Point
# ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Quick-start: just run with defaults
    #   $ python nhanes_stress_mortality.py
    #
    # Custom config:
    #   cfg = NHANESConfig(
    #       cycles=["2015-2016", "2017-2018"],  # fewer cycles = faster
    #       quantile_cutoff=0.85,                # top 15%
    #   )
    #   pipeline = AllostaticLoadPipeline(cfg).run()

    pipeline = AllostaticLoadPipeline().run()
