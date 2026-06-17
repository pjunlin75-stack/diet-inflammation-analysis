# 膳食因素与全身炎症分析
### NHANES 2017-2018 数据分析项目 (Python)

![Omega Ratio vs CRP](results/omega_crp_result.png)

这是一个**可复用模板仓库**，用于探索不同膳食摄入与 hs-CRP（全身炎症标志物）的关系。

## 项目目标
- 学习 Python 数据分析流程（下载、清洗、统计、可视化）
- 探索营养流行病学经典问题
- 建立个人 GitHub Portfolio

## 已完成分析

| 脚本 | 变量 | 主要发现 | 文件 |
|------|------|----------|------|
| `omega_crp_nhanes.py` | Omega-6:Omega-3 比例 | 无显著关联 | [脚本](scripts/omega_crp_nhanes.py) |
| `sodium_crp_analysis.py` | 膳食钠摄入 | 弱负相关 (p<0.001) | [脚本](scripts/sodium_crp_analysis.py) |
| `sugar_crp_analysis.py` | 总糖摄入 | 无显著关联 | [脚本](scripts/sugar_crp_analysis.py) |
| `fiber_crp_analysis.py` | 膳食纤维 (可选) | - | - |
| `metabolic_depression_nhanes.py` | 代谢综合征 → 抑郁 | OR=1.68 (p<0.0001) | [脚本](scripts/metabolic_depression_nhanes.py) |

## 如何使用模板
1. 复制 `scripts/` 里的任意 `.py` 文件
2. 修改 `variable_code` 和标题
3. 运行即可得到新分析

## 运行环境
```bash
cd [你的文件夹]
source myenv/bin/activate
python scripts/omega_crp_nhanes.py
