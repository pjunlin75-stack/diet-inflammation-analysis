"""
sodium_crp_analysis.py - 优化版
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'Hei', 'Songti SC']
plt.rcParams['axes.unicode_minus'] = False

print("=== 膳食钠摄入与 CRP 分析（优化版）===\n")

demo = pd.read_sas("DEMO_J.xpt")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
diet = pd.read_sas("DR1TOT_J.xpt")[["SEQN", "DR1TSODI"]]
crp = pd.read_sas("HSCRP_J.xpt")[["SEQN", "LBXHSCRP"]]

df = demo.merge(diet, on="SEQN").merge(crp, on="SEQN")

df = df[
    (df.RIDAGEYR >= 18) &
    (df.DR1TSODI > 0) & 
    (df.LBXHSCRP.notna())
].copy()

df["logCRP"] = np.log(df.LBXHSCRP)
df["female"] = (df.RIAGENDR == 2).astype(int)

print(f"最终分析样本: n = {len(df)}")
print(f"  中位钠摄入 (mg/天) = {df.DR1TSODI.median():.0f}")

r, p = spearmanr(df.DR1TSODI, df.LBXHSCRP)
print(f"\nSpearman 相关: r = {r:.3f}, p = {p:.4f}")

df["quintile"] = pd.qcut(df.DR1TSODI, 5, labels=["Q1(最低)", "Q2", "Q3", "Q4", "Q5(最高)"])

summary = df.groupby("quintile", observed=True).agg(
    median_CRP=("LBXHSCRP", "median"),
    median_sodium=("DR1TSODI", "median"),
    n=("SEQN", "size")
)
print("\n钠摄入分位数 vs 中位 CRP:")
print(summary.round(2))

# 画图
fig, ax = plt.subplots(figsize=(9, 6))
colors = ['#2ecc71', '#82e0aa', '#f9e79f', '#f5b041', '#e74c3c']

bars = ax.bar(range(5), summary.median_CRP.values, color=colors, edgecolor='black')

ax.set_xticks(range(5))
ax.set_xticklabels(summary.index)
ax.set_xlabel("膳食钠摄入量（五分位）", fontsize=12)
ax.set_ylabel("中位 hs-CRP (mg/L)", fontsize=12)
ax.set_title("膳食钠摄入与全身炎症 (hs-CRP)\nNHANES 2017-2018", fontsize=14, pad=20)

# 添加数值标签
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
            f'{height:.2f}', ha='center', va='bottom', fontsize=11)

plt.savefig("sodium_crp_result.png", dpi=200, bbox_inches='tight')
print("\n✅ 图表已保存: sodium_crp_result.png")
print("   用 Finder 打开 Downloads 文件夹即可查看")
