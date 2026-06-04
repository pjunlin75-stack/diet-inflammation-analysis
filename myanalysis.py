"""
my_analysis.py
研究问题：膳食【XXX】与 CRP 的关系
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

# ================== 只改这里 1、2、3 处 ==================

# 【1】改变量名和中文名称
variable_code = "DR1TSUGR"          # ← 改成你要的变量代码（如 DR1TFIBE）
variable_name = "总糖摄入 (g/天)"    # ← 改成中文描述

# 【2】改图表标题
title = "膳食总糖摄入与全身炎症 (hs-CRP)"

# 【3】（可选）改文件名
output_image = "sugar_crp_result.png"

# =========================================================

print(f"=== {title} 分析 ===\n")

demo = pd.read_sas("DEMO_J.xpt")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
diet = pd.read_sas("DR1TOT_J.xpt")[["SEQN", variable_code]]
crp = pd.read_sas("HSCRP_J.xpt")[["SEQN", "LBXHSCRP"]]

df = demo.merge(diet, on="SEQN").merge(crp, on="SEQN")

df = df[
    (df.RIDAGEYR >= 18) &
    (df[variable_code] > 0) & 
    (df.LBXHSCRP.notna())
].copy()

print(f"最终分析样本: n = {len(df)}")
print(f"  中位{variable_name} = {df[variable_code].median():.2f}")

r, p = spearmanr(df[variable_code], df.LBXHSCRP)
print(f"\nSpearman 相关: r = {r:.3f}, p = {p:.4f}")

df["quintile"] = pd.qcut(df[variable_code], 5, 
                        labels=["Q1(最低)", "Q2", "Q3", "Q4", "Q5(最高)"])

summary = df.groupby("quintile", observed=True).agg(
    median_CRP=("LBXHSCRP", "median"),
    median_var=(variable_code, "median"),
    n=("SEQN", "size")
)
print(f"\n{variable_name}分位数 vs 中位 CRP:")
print(summary.round(2))

# 画图
fig, ax = plt.subplots(figsize=(9, 6))
colors = ['#2ecc71', '#82e0aa', '#f9e79f', '#f5b041', '#e74c3c']

bars = ax.bar(range(5), summary.median_CRP.values, color=colors)

ax.set_xticks(range(5))
ax.set_xticklabels(summary.index)
ax.set_xlabel(variable_name, fontsize=12)
ax.set_ylabel("中位 hs-CRP (mg/L)", fontsize=12)
ax.set_title(f"{title}\nNHANES 2017-2018", fontsize=14, pad=20)

plt.savefig(output_image, dpi=200, bbox_inches='tight')
print(f"\n✅ 图表已保存: {output_image}")
