# 膳食因素与全身炎症分析模板
### NHANES 2017-2018 (Python)

![Omega-6:3 vs CRP](results/omega_crp_result.png)

这是一个**可复用**的营养流行病学分析模板仓库，用于探索不同膳食摄入与 hs-CRP（炎症标志物）的关系。

## 项目亮点

- 从零基础开始，用 AI 辅助完成
- 包含多个真实变量分析（Omega、钠、糖等）
- 代码结构清晰，适合作为学习模板

## 文件结构
├── scripts/                  # 分析脚本
│   ├── omega_crp_nhanes.py
│   ├── sodium_crp_analysis.py
│   ├── sugar_crp_analysis.py
│   └── fiber_crp_analysis.py
├── results/                  # 生成的图表
│   ├── *.png
└── README.md
text## 如何使用本模板

1. 复制任意 `scripts/` 里的 `.py` 文件
2. 修改 `variable_code` 和标题
3. 运行即可得到新分析结果

## 主要发现（示例）

- **Omega-6:Omega-3 比例**：无显著关联
- **膳食钠**：呈现弱负相关（p<0.001）
- **总糖摄入**：无显著关联

## 技术栈

- Python 3
- pandas, numpy, scipy, matplotlib
- NHANES 公开数据

---
