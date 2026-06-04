"""
sugar_crp_analysis.py
研究问题：膳食总糖摄入与全身炎症（hs-CRP）的关系
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

print("=== 膳食总糖摄入与 CRP 分析 ===\n")

# 1. 读取数据
demo = pd.read_sas("DEMO_J.xpt")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
diet = pd.read_sas("DR1TOT_J.xpt")[["SEQN", "DR1TSUGR"]]   # ← 总糖 DR1TSUGR
crp = pd.read_sas("HSCRP_J.xpt")[["SEQN", "LBXHSCRP"]]

# 2. 合并数据
df = demo.merge(diet, on="SEQN").merge(crp, on="SEQN")

# 3. 数据清洗
df = df[
    (df.RIDAGEYR >= 18) &
    (df.DR1TSUGR > 0) & 
    (df.LBXHSCRP.notna())
].copy()

df["logCRP"] = np.log(df.LBXHSCRP)
df["female"] = (df.RIAGENDR == 2).astype(int)

print(f"最终分析样本: n = {len(df)}")
print(f"  中位总糖摄入 (g/天) = {df.DR1TSUGR.median():.1f}")

# 4. Spearman 相关
r, p = spearmanr(df.DR1TSUGR, df.LBXHSCRP)
print(f"\nSpearman 相关: r = {r:.3f}, p = {p:.4f}")

# 5. 分位数分析
df["quintile"] = pd.qcut(df.DR1TSUGR, 5, labels=["Q1(最低)", "Q2", "Q3", "Q4", "Q5(最高)"])

summary = df.groupby("quintile", observed=True).agg(
    median_CRP=("LBXHSCRP", "median"),
    median_sugar=("DR1TSUGR", "median"),
    n=("SEQN", "size")
)
print("\n总糖摄入分位数 vs 中位 CRP:")
print(summary.round(2))

# 6. 画图
fig, ax = plt.subplots(figsize=(9, 6))
colors = ['#2ecc71', '#82e0aa', '#f9e79f', '#f5b041', '#e74c3c']

bars = ax.bar(range(5), summary.median_CRP.values, color=colors, edgecolor='black')

ax.set_xticks(range(5))
ax.set_xticklabels(summary.index)
ax.set_xlabel("膳食总糖摄入量（五分位）", fontsize=12)
ax.set_ylabel("中位 hs-CRP (mg/L)", fontsize=12)
ax.set_title("膳食总糖摄入与全身炎症 (hs-CRP)\nNHANES 2017-2018", fontsize=14, pad=20)

# 添加数值标签
for i, bar in enumerate(bars):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
            f'{height:.2f}', ha='center', va='bottom', fontsize=11)

plt.savefig("sugar_crp_result.png", dpi=200, bbox_inches='tight')
print("\n✅ 图表已保存: sugar_crp_result.png")
print("   建议用命令打开: open sugar_crp_result.png")
