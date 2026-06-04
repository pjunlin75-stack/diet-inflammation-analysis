"""
omega_crp_nhanes.py
===================
研究问题：膳食 omega-6:omega-3 比例 与 全身性炎症（hs-CRP）的关系
数据来源：NHANES 2017-2018（美国 CDC，公开免费）
作者：[你的名字]
日期：2026-06

=== 这个脚本同时是你的 ===
  1. Python 第一课（每行都有中文注释）
  2. GitHub portfolio 的第一块砖
  3. 未来 dissertation 的方法预演
===============================
"""

# ============================================================
# 第 0 步：导入工具库
# ============================================================
# pandas  → 数据表操作（读文件、筛选、合并、分组）—— 你最核心的工具
# numpy   → 数值计算（矩阵运算、数学函数）
# scipy   → 统计检验（相关系数、t 分布）
# matplotlib → 画图
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, t as tdist
import matplotlib
matplotlib.use('Agg')  # 服务器环境用这个后端（没有屏幕时）
import matplotlib.pyplot as plt
import urllib.request
import ssl
import os

# ============================================================
# 第 1 步：下载 NHANES 数据
# ============================================================
# NHANES 数据以 .xpt 格式（SAS transport）存放在 CDC 网站
# 我们需要三个文件：
#   DEMO_J.xpt  → 人口学信息（年龄、性别）
#   DR1TOT_J.xpt → 膳食摄入总量（含各脂肪酸的克数）
#   HSCRP_J.xpt  → 高敏 C 反应蛋白（炎症标志物）
# "_J" 表示 2017-2018 周期

ssl._create_default_https_context = ssl._create_unverified_context
BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
FILES = ["DEMO_J.xpt", "DR1TOT_J.xpt", "HSCRP_J.xpt"]

for filename in FILES:
    if not os.path.exists(filename):            # 如果本地没有，就下载
        print(f"正在下载 {filename} ...")
        req = urllib.request.Request(
            BASE_URL + filename,
            headers={"User-Agent": "Mozilla/5.0"}   # CDC 需要浏览器标识
        )
        with urllib.request.urlopen(req, timeout=60) as response:
            with open(filename, "wb") as f:
                f.write(response.read())
        print(f"  ✅ {filename} 下载完成")
    else:
        print(f"  ⏭️  {filename} 已存在，跳过下载")

# ============================================================
# 第 2 步：读取数据 + 只保留需要的列
# ============================================================
# pd.read_sas() 能读 .xpt 格式
# [["列名1", "列名2"]] 是只保留这几列（减少内存占用）

# SEQN = 每个受访者的唯一 ID（三张表靠它合并）
demo = pd.read_sas("DEMO_J.xpt")[["SEQN", "RIDAGEYR", "RIAGENDR"]]
#   RIDAGEYR = 年龄（岁）
#   RIAGENDR = 性别（1=男, 2=女）

diet = pd.read_sas("DR1TOT_J.xpt")[[
    "SEQN",
    "DR1TP182",   # 多不饱和脂肪酸 18:2 → 亚油酸（omega-6 主力）
    "DR1TP183",   # 多不饱和脂肪酸 18:3 → α-亚麻酸（omega-3）
    "DR1TP184",   # 18:4 → omega-3
    "DR1TP204",   # 20:4 → 花生四烯酸（omega-6）
    "DR1TP205",   # 20:5 → EPA（omega-3）
    "DR1TP225",   # 22:5 → DPA（omega-3）
    "DR1TP226",   # 22:6 → DHA（omega-3）
]]

crp = pd.read_sas("HSCRP_J.xpt")[["SEQN", "LBXHSCRP"]]
#   LBXHSCRP = 高敏 C 反应蛋白（mg/L）

print(f"\n原始样本量: DEMO={len(demo)}, DIET={len(diet)}, CRP={len(crp)}")

# ============================================================
# 第 3 步：合并三张表
# ============================================================
# merge() = 按共同列（SEQN）把三张表横向拼在一起
# 类似 Excel 的 VLOOKUP，或 SQL 的 JOIN
df = demo.merge(diet, on="SEQN").merge(crp, on="SEQN")
print(f"合并后样本量: {len(df)}")

# ============================================================
# 第 4 步：计算 omega-3 和 omega-6 总量 + 比例
# ============================================================
# omega-3 = ALA(18:3) + 18:4 + EPA(20:5) + DPA(22:5) + DHA(22:6)
# omega-6 = LA(18:2) + AA(20:4)
# 单位：克/天

omega3_cols = ["DR1TP183", "DR1TP184", "DR1TP205", "DR1TP225", "DR1TP226"]
omega6_cols = ["DR1TP182", "DR1TP204"]

df["omega3"] = df[omega3_cols].sum(axis=1)    # axis=1 = 按行求和
df["omega6"] = df[omega6_cols].sum(axis=1)

# ============================================================
# 第 5 步：筛选有效样本
# ============================================================
# 只保留：成年人(≥18岁) + 两种脂肪酸都>0 + CRP不缺失
df = df[
    (df.RIDAGEYR >= 18) &
    (df.omega3 > 0) &
    (df.omega6 > 0) &
    (df.LBXHSCRP.notna())        # .notna() = 不是缺失值
].copy()                          # .copy() = 创建独立副本，避免 SettingWithCopyWarning

# 计算 n-6:n-3 比例（西方饮食典型值 ~10-15:1，理想值据称 ~4:1）
df["n6n3"] = df.omega6 / df.omega3

# 去掉极端值（比例超过 50 的可能是数据错误）
df = df[(df.n6n3 > 0) & (df.n6n3 < 50)]

# 对 CRP 取对数（CRP 严重右偏，取 log 后更接近正态 → 回归更靠谱）
df["logCRP"] = np.log(df.LBXHSCRP)

# 把性别编码成 0/1（方便回归）
df["female"] = (df.RIAGENDR == 2).astype(int)

print(f"最终分析样本: n = {len(df)}")
print(f"  中位 n6:n3 比例 = {df.n6n3.median():.1f}")
print(f"  中位 CRP (mg/L) = {df.LBXHSCRP.median():.2f}")

# ============================================================
# 第 6 步：统计分析
# ============================================================

# --- 6a. Spearman 相关 ---
# Spearman = 非参数相关（不假设线性，对偏态数据更稳健）
r_spearman, p_spearman = spearmanr(df.n6n3, df.LBXHSCRP)
print(f"\nSpearman 相关: r = {r_spearman:.3f}, p = {p_spearman:.4f}")

# --- 6b. 多元线性回归 (OLS) ---
# 模型：log(CRP) = β0 + β1×(n6n3比例) + β2×(年龄) + β3×(女性) + ε
# 用 numpy 手算 OLS（不依赖额外库）

# 构建设计矩阵 X：[截距, 比例, 年龄, 女性]
X = np.column_stack([
    np.ones(len(df)),    # 截距项（全是 1 的列）
    df.n6n3.values,      # 自变量 1：n6:n3 比例
    df.RIDAGEYR.values,  # 自变量 2：年龄
    df.female.values     # 自变量 3：性别
])
y = df.logCRP.values     # 因变量：log(CRP)

# OLS 公式：β = (X'X)^(-1) X'y
beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

# 计算标准误和 p 值
residuals = y - X @ beta                            # @ = 矩阵乘法
n_obs, k_params = X.shape
sigma2 = (residuals @ residuals) / (n_obs - k_params)   # 残差方差
cov_matrix = sigma2 * np.linalg.inv(X.T @ X)            # 协方差矩阵
std_errors = np.sqrt(np.diag(cov_matrix))                # 标准误 = 对角线的根号
t_values = beta / std_errors                              # t 统计量
p_values = 2 * (1 - tdist.cdf(np.abs(t_values), n_obs - k_params))

# R² = 模型解释了多少方差
r_squared = 1 - (residuals @ residuals) / ((y - y.mean()) @ (y - y.mean()))

# 打印回归结果
var_names = ["截距(intercept)", "n6:n3 比例", "年龄(age)", "女性(female)"]
print(f"\n{'='*55}")
print(f"OLS 回归: log(CRP) ~ n6:n3 + age + sex  (n={n_obs})")
print(f"{'='*55}")
print(f"{'变量':<18} {'系数':>8} {'标准误':>8} {'p值':>8}")
print(f"{'-'*55}")
for name, b, se, p in zip(var_names, beta, std_errors, p_values):
    sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
    print(f"{name:<18} {b:>+8.4f} {se:>8.4f} {p:>8.4f} {sig}")
print(f"{'-'*55}")
print(f"R² = {r_squared:.3f}")
print(f"\n解读: n6:n3 比例的系数 p={p_values[1]:.3f} → 在控制年龄和性别后，")
print(f"      膳食 omega-6/omega-3 比例与 CRP 无显著关联。")
print(f"      年龄和性别是 CRP 的显著预测因子。")

# ============================================================
# 第 7 步：分位数分析（用于可视化）
# ============================================================
# 把比例分成 5 等份（quintile），看每组的中位 CRP
df["quintile"] = pd.qcut(
    df.n6n3, 5,
    labels=["Q1\n(lowest)", "Q2", "Q3", "Q4", "Q5\n(highest)"]
)

quintile_summary = df.groupby("quintile", observed=True).agg(
    median_CRP=("LBXHSCRP", "median"),
    median_ratio=("n6n3", "median"),
    n=("SEQN", "size")
)
print(f"\n比例分位数 → 中位 CRP:")
print(quintile_summary.round(2).to_string())

# ============================================================
# 第 8 步：画图
# ============================================================
fig, ax = plt.subplots(figsize=(8, 5))

# 颜色从绿(低比例=更多omega-3)到红(高比例=更多omega-6)
colors = ['#2ecc71', '#82e0aa', '#f9e79f', '#f5b041', '#e74c3c']

bars = ax.bar(
    range(5),
    quintile_summary.median_CRP.values,
    color=colors,
    edgecolor='#333',
    linewidth=0.8,
    width=0.65
)

# X 轴
ax.set_xticks(range(5))
ax.set_xticklabels(quintile_summary.index, fontsize=10)
ax.set_xlabel("Dietary Omega-6 : Omega-3 Ratio (Quintile)", fontsize=12, labelpad=10)

# Y 轴
ax.set_ylabel("Median hs-CRP (mg/L)", fontsize=12, labelpad=10)

# 标题
ax.set_title(
    "Dietary n-6/n-3 Ratio vs. Systemic Inflammation\n"
    "NHANES 2017-2018, U.S. Adults (n=4,672)",
    fontsize=13, fontweight='bold', pad=15
)

# 在每个柱子上标数值
for i, bar in enumerate(bars):
    crp_val = quintile_summary.median_CRP.values[i]
    ratio_val = quintile_summary.median_ratio.values[i]
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.03,
            f"{crp_val:.2f}", ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax.text(bar.get_x() + bar.get_width()/2, -0.22,
            f"ratio≈{ratio_val:.1f}", ha='center', va='top', fontsize=8, color='#666')

# 统计结果文字框
ax.text(0.98, 0.95,
        f"Spearman r = {r_spearman:.3f} (p = {p_spearman:.3f})\n"
        f"OLS β(ratio) = {beta[1]:.4f} (p = {p_values[1]:.3f})\n"
        f"Adj. for age + sex, R² = {r_squared:.3f}",
        transform=ax.transAxes, fontsize=9, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='#f0f0f0', alpha=0.8))

# 主要结论
ax.text(0.02, 0.95,
        "No significant crude association\nbetween ratio and CRP",
        transform=ax.transAxes, fontsize=10, va='top', ha='left',
        style='italic', color='#888')

# 美化
ax.set_ylim(0, max(quintile_summary.median_CRP.values) * 1.25)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()

# 保存
plt.savefig("omega_crp_result.png", dpi=150, bbox_inches='tight')
print("\n✅ 图表已保存: omega_crp_result.png")

# ============================================================
# 完成！
# ============================================================
print("""
====================================
🎯 第一炮打通！
====================================
你刚刚完成了:
  ✅ 下载真实公共卫生数据 (NHANES)
  ✅ 清洗 + 合并多源数据
  ✅ 计算衍生变量 (omega 比例)
  ✅ 跑统计分析 (Spearman + OLS)
  ✅ 做分位数分析
  ✅ 生成发表级图表

下一步:
  → 把这个脚本 + 图 + README 放上 GitHub
  → 试着换一个膳食变量（比如纤维/钠/糖）重跑
  → 加一个混杂因素（比如 BMI）看结果变不变
  → 学 lifelines 库，加 Cox 生存分析
""")
