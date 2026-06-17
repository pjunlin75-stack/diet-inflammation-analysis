# Chronic Stress and Mortality Risk – NHANES Survival Analysis
### Personal Educational Project (Python)

This repository explores the relationship between chronic psychological stress (measured by PHQ-9) and long-term all-cause mortality using publicly available NHANES data (2005–2018).

**Disclaimer**: This is a **personal learning project** for educational and portfolio purposes only. It is not peer-reviewed research, not medical advice, and should not be used for clinical decisions.

## Project Goals
- Practice end-to-end data analysis pipeline (download, cleaning, feature engineering, survival modeling)
- Understand how psychological stress may relate to physiological wear-and-tear and mortality
- Build reusable Python code for public health microdata analysis

## Key Components
- Automated NHANES data loader (multiple cycles)
- PHQ-9 scoring and high-stress phenotype definition
- Metabolic/inflammatory biomarkers (CRP, blood pressure)
- Cox Proportional Hazards survival model
- Kaplan-Meier survival curves and dose-response visualization

## Repository Structure
```
├── scripts/           # Analysis scripts
├── results/           # Generated plots
├── output/            # Analytic datasets (optional)
└── README.md
```

## How to Run
```bash
cd [project folder]
source myenv/bin/activate
python nhanes_stress_mortality.py
```

## Technologies Used
- Python 3
- pandas, numpy, lifelines (survival analysis), matplotlib
- Public NHANES + Linked Mortality data (CDC)

## Learning Outcomes
- From zero Python background to building a full survival analysis pipeline
- Practical experience with complex survey data and time-to-event modeling
```

