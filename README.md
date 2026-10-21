# scRNA-seq 自动化分析流程 (GEO单细胞数据处理)

[![R](https://img.shields.io/badge/R-%3E%3D%204.1.0-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v4%20%2F%20v5-green.svg)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

本仓库基于 **R** 语言开发，是一套专门针对 **GEO (Gene Expression Omnibus)** 数据库下载的单细胞转录组（scRNA-seq）数据设计的自动化分析流程。

公共数据库（尤其是 GEO）的单细胞数据格式五花八门（如 10x 标准的三件套、单文件 `.txt`/`.csv` 矩阵、`.tar` 或 `.gz` 压缩包等），手动读取和清洗极其耗时。本仓库提供**一键式处理脚本**，只需将下载的数据放入指定文件夹，即可自动完成从数据读取、质控、降维聚类到差异分析与功能富集的完整常规流程。

---

## ✨ 核心特色

- 📦 **GEO 数据无缝对接**：内置自动解压机制，兼容多种主流单细胞矩阵格式（10x Genomics 标准格式、单文件表达矩阵等）。
- 🚀 **常规数据一键处理**：提供主控脚本，自动执行质控、整合（Harmony）、聚类、注释及可视化，直接输出可用图表。
- 🧩 **高度模块化**：底层拆分为独立的 R 脚本，既支持全自动跑完，也方便随时暂停，对特定参数（如 `resolution`、`dims`）进行微调。
- 🎯 **面向机制分析**：内置基因集评分（如衰老、SASP、炎症、代谢）和差异分析模块，方便直接生成用于撰写文章的核心图表。

---

## 📁 仓库结构与模块说明

### 📥 1. `scRNA_sort` — GEO 数据读取与预处理
- **自动解压**：识别并解压 GEO 下载的 `.tar`, `.gz` 等压缩包。
- **格式兼容**：支持 10x 标准格式（`barcodes.tsv`, `features.tsv`, `matrix.mtx`）及常见的单文件 mRNA 表达矩阵。
- **数据整理**：自动处理样本文件名称，处理后支持scRNA_QC直接输入。

### 🛠️ 2. `scRNA_QC` — 一键质控
- 自动计算线粒体基因比例（`percent.mt`）。
- 绘制质控小提琴图（Violin Plot），支持自动调整过滤阈值。

### 📊 3. `scRNA_Scale` — 降维与聚类
- 数据标准化（`NormalizeData`、 `FindVariableFeatures`或 `ScaleData`）。
- 自动提取高变基因（HVGs），进行 PCA 降维、 UMAP 降维和 t-SNE 降维。
- 内置 **Harmony** 算法，若包含多个 GEO 样本，可一键去除批次效应。
- 使用 SingleR 和 ScType 进行初步的自动注释。

### 🧬 4. `scRNA_DEG` — 差异表达分析
包含两个子模块：
- **`Celltype_DEG`**：用于寻找各细胞亚群的特异性 Marker 基因（支持一键输出 DotPlot 和 FeaturePlot）。
- **`Group_DEG`**：用于比较实验组与对照组（如 OA vs. Control）之间的差异基因。

### 📈 5. `scRNA_enrich` — GO / KEGG 功能富集
- 针对差异基因，自动进行 GO（BP/CC/MF）、 KEGG 和 GSEA 富集。
- 生成条形图（Barplot）和气泡图（Dotplot），并输出富集结果表格。

### 🎯 6. `Set_score` — 通路/基因集打分
- 基于 `AUCell` 或 `AddModuleScore`，对每个细胞计算特定功能状态（如**衰老、SASP状态、炎症、代谢**）的活性评分。
- 直观展示不同细胞亚群或不同分组间的功能状态差异。

### 🔗 7. 其他辅助模块
- **`scRNA_cor`**：计算基因与基因、基因与通路评分之间的相关性。

---

## 🚀 快速开始 (Quick Start)

### 第一步：准备数据
将 GEO 下载的原始数据文件夹（例如 `GSE123456_RAW.tar`）放入你的工作目录。
若为10x 标准格式（`barcodes.tsv`, `features.tsv`, `matrix.mtx`）格式，使用 scRNA_cor 进行处理，其他格式需自行处理。

### 第二步：
生成样本信息文件，要求第一列为样本，第二列为样本所属分组，第三列为差异比较组，对照__处理。

```r
# 基础配置示例
data_dir <- "path/to/your/GEO_data"   # GEO 数据存放路径
species <- "human"                    # "human" 或 "mouse"
project_name <- "My_GEO_Project"      # 项目名称（用于输出文件命名）
