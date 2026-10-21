# scRNA-seq 自动化分析流程 (GEO 单细胞数据处理)

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

本仓库基于 **R** 语言开发，是一套专门针对 **GEO (Gene Expression Omnibus)** 数据库下载的单细胞转录组（scRNA-seq）数据设计的自动化分析流程。

公共数据库（尤其是 GEO）的单细胞数据格式五花八门（如 10x 标准的三件套、单文件 `.txt`/`.csv` 矩阵、`.tar` 或 `.gz` 压缩包等），手动读取和清洗极其耗时。本仓库提供**一键式处理脚本**，只需将下载的数据放入指定文件夹，即可自动完成从数据读取、质控、降维聚类、**细胞类型注释**到差异分析与功能富集的完整常规流程。

## ✨ 核心特色

- **GEO 数据无缝对接**：内置自动解压机制，兼容多种主流单细胞矩阵格式（10x Genomics 标准格式、单文件表达矩阵等）。
- **常规数据一键处理**：提供主控脚本，自动执行质控、整合（Harmony）、聚类、**细胞注释**及可视化，直接输出可用图表。
- **高度模块化**：底层拆分为独立的 R 脚本，既支持全自动跑完，也方便随时暂停，对特定参数（如 `resolution`、`dims`）进行微调。
- **面向机制分析**：内置基因集评分（如衰老、SASP、炎症、代谢）和差异分析模块，方便直接生成用于撰写文章的核心图表。

## 📁 仓库结构与模块说明

### 1. `scRNA_sort` — GEO 数据读取与预处理
- **自动解压**：识别并解压 GEO 下载的 `.tar`, `.gz` 等压缩包。
- **格式兼容**：支持 10x 标准格式（`barcodes.tsv`, `features.tsv`, `matrix.mtx`）及常见的单文件 mRNA 表达矩阵。
- **数据整理**：自动处理样本文件名称，处理后支持 `scRNA_QC` 直接输入。

### 2. `scRNA_QC` — 一键质控
- 自动计算线粒体基因比例（`percent.mt`）。
- 绘制质控小提琴图（Violin Plot），支持自动调整过滤阈值。

### 3. `scRNA_Scale` — 降维与聚类
- 数据标准化（`NormalizeData`、`FindVariableFeatures` 和 `ScaleData`）。
- 自动提取高变基因（HVGs），进行 PCA 降维、UMAP 降维和 t-SNE 降维。
- 内置 **Harmony** 算法，若包含多个 GEO 样本，可一键去除批次效应。
- 使用 **SingleR** 和 **ScType** 进行初步注释，为后续基于 Marker 基因注释提供参考。

### 4. `CellType_annotate` — 细胞类型注释
针对细胞注释环节，使用 **SingleR**、**ScType** 并结合 Marker 基因的策略，确保注释结果的准确性与可重复性：
- **基于 Marker 基因的手动注释**：通过 `FindAllMarkers()` 鉴定各细胞簇的特异性高表达基因，结合 `DotPlot` 和 `FeaturePlot` 可视化验证，再与已知经典 Marker 基因（如 `CD3E`、`PTPRC`、`COL2A1` 等）比对，完成最终的细胞类型指派。

> **建议流程**：先用 SingleR + ScType 获得初步注释，再通过手动 Marker 基因验证和修正，三者交叉验证以提高注释可靠性。

### 5. `scRNA_DEG` — 差异表达分析
包含两个子模块：
- **`Celltype_DEG`**：用于寻找各细胞亚群的特异性 Marker 基因（支持一键输出 DotPlot 和 FeaturePlot）。
- **`Group_DEG`**：用于比较实验组与对照组（如 OA vs. Control）之间的差异基因。

### 6. `scRNA_enrich` — GO / KEGG 功能富集
- 针对差异基因，自动进行 GO（BP/CC/MF）、KEGG 和 GSEA 富集。
- 生成条形图（Barplot）和气泡图（Dotplot），并输出富集结果表格。

### 7. `Set_score` — 通路/基因集打分
- 基于 `AUCell` 或 `AddModuleScore`，对每个细胞计算特定功能状态（如**衰老、SASP状态、炎症、代谢**）的活性评分。
- 直观展示不同细胞亚群或不同分组间的功能状态差异。

### 8. 其他辅助模块
- **`scRNA_cor`**：计算基因与基因、基因与通路评分之间的相关性。
- **`Gene_plot`**：针对特定选出的基因，进行组间和细胞类型间的多水平可视化。

## 🚀 快速开始 (Quick Start)

### 第一步：准备数据
将 GEO 下载的原始数据文件夹（例如 `GSE123456_RAW.tar`）放入你的工作目录。若为 10x 标准格式（`barcodes.tsv`, `features.tsv`, `matrix.mtx`），使用 `scRNA_sort` 进行处理；其他格式需自行整理。

### 第二步：生成样本信息文件
要求第一列为样本名，第二列为样本所属分组，第三列为差异比较组（对照组在前，如 `Control__OA`）。

### 第三步：查看帮助文档并运行
依据脚本帮助文档，参考脚本运行所需的必须参数与可选参数，补全参数后进行分析。

```bash
scRNA_Scale --help
```

## 💻 配置环境

建议使用 R (>= 4.1.0) 环境，核心依赖包如下：

```r
# 核心处理与可视化
install.packages(c("Seurat", "tidyverse", "harmony", "ggplot2", "patchwork"))

# 细胞注释
install.packages("SingleR")
BiocManager::install("celldex")        # SingleR 参考数据集

# 基因集评分与富集分析
BiocManager::install(c("AUCell", "clusterProfiler", "org.Hs.eg.db"))
# 注意：如果数据为小鼠，请将 org.Hs.eg.db 替换为 org.Mm.eg.db
```

> `ScType` 的安装请参考其 GitHub 仓库（`IanevskiAleksandr/sc-type`）。

## ⚠️ 注意事项

- **物种注释**：GEO 数据如果未标明物种，请通过基因命名特征（如人类全大写 `CD3E`，小鼠首字母大写 `Cd3e`）确认，并在 `scRNA_enrich` 中设置准确的数据库。
- **批次效应**：如果下载的 GEO 数据由不同平台或多个独立实验组成，强烈建议在 `scRNA_Scale` 阶段开启 `Harmony` 降维，并在 UMAP 图中检查不同样本是否混合均匀。
- **细胞注释的交叉验证**：`SingleR` 和 `ScType` 的结果是“参考意见”而非“最终答案”。务必通过 `FindAllMarkers` 输出的 Marker 基因进行手动校验，尤其注意检查自动注释结果与经典 Marker 基因表达是否一致。对于自动注释置信度较低的细胞簇，应结合 `DotPlot` 和 `FeaturePlot` 逐簇确认。
- **聚类分辨率**：若注释后发现某些细胞簇无法被明确归类，或 Marker 基因表达模式混杂，可尝试调整 `FindClusters` 的 `resolution` 参数重新聚类。

## 📄 许可协议

本仓库脚本仅供学习与科研参考，请遵循相关数据的使用规范与开源协议。
