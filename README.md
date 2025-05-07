# Amplicaon_Analysis
## 扩增子分析
### Amplicaon_processing.py 脚本功能与用途

该脚本用于自动化处理扩增子（amplicon）测序数据，从原始双端 FASTQ 数据到最终 OTU/ASV 表，并可选地进行分类注释。主要流程包括：  
1. 样本解复用（demultiplex）  
2. 双端合并（merge paired‐ends）  
3. 质量过滤（quality filtering）  
4. 序列去重复（dereplication）  
5. OTU/ASV 聚类（clustering）  
6. 嵌合体检测（chimera detection）  
7. OTU 表生成（OTU table）  
8. 可选的 SINTAX 分类注释（taxonomy assignment）

---

## 输入

### 一、命令行参数  
| 参数                | 说明                                                           | 示例                               |
|---------------------|----------------------------------------------------------------|------------------------------------|
| `-i, --input_dir`   | 包含所有原始测序文件（FASTQ）的目录（绝对或相对路径）。       | `raw_data/`                        |
| `-m, --metadata_file` | 元数据文件路径（TSV 格式，无表头）。                          | `metadata.tsv`                     |
| `-o, --output_dir`  | 结果输出目录（若不存在则自动创建）。                            | `results/`                         |
| `-t, --threads`     | 并行使用的线程数（整数）。                                     | `4`                                | 

### 二、元数据文件格式（TSV，**无表头**）  
脚本假设元数据文件按以下列顺序排列：  
- **run_id**：测序运行编号  
- **sample_id**：样本编号，用于命名中间及输出文件  
- **forward_primer**：正向引物序列（用于解复用）  
- **reverse_primer**：反向引物序列（用于解复用）  
- **forward_file**：原始正向 FASTQ 文件名（相对于 `input_dir`）  
- **reverse_file**：原始反向 FASTQ 文件名（相对于 `input_dir`）  

#### 示例内容  
```tsv
run01    SampleA    AGAGTTTGATCMTGGCTCAG    CGGTTACCTTGTTACGACTT    SampleA_R1.fastq.gz    SampleA_R2.fastq.gz
run01    SampleB    AGAGTTTGATCMTGGCTCAG    CGGTTACCTTGTTACGACTT    SampleB_R1.fastq.gz    SampleB_R2.fastq.gz
```

### 三、交互式参数  
脚本在运行过程中会根据用户选择提示：  
1. **PCR 产物长度类型**（固定范围或指定多个值）及相应的长度数值  
2. **聚类策略**：选择 UPARSE3（OTU）或 UNOISE3（ASV）  
3. **嵌合体检测方法**：de novo 或 参考数据库  
4. **是否执行 SINTAX 分类注释**（yes/no）及分类数据库路径（若选择 yes）

---

## 输出

脚本在指定的 `output_dir` 下，按照每个处理步骤创建子目录并生成相应文件：

| 步骤                  | 子目录               | 主要输出文件示例                           |
|-----------------------|----------------------|--------------------------------------------|
| 1. 解复用             | `1-demultiplex/`     | `SampleA.R1.fastq`, `SampleA.R2.fastq`     |
| 2. 合并               | `2-merge/`           | `SampleA.merged.fastq`                     |
| 3. 质量过滤           | `3-quality/`         | `SampleA.filtered.fasta`                   |
| 4. 去重复             | `4-dereplicate/`     | `SampleA.derep.fasta`, `all_samples_derep.fasta` |
| 5. 聚类               | `5-cluster/`         | `otus.fasta` (或 centroids.fasta)          |
| 6. 嵌合体检测         | `6-chimera/`         | `otus_nochim.fasta`                        |
| 7. OTU 表生成         | `7-OTU/`             | `otu_table.txt`                            |
| 8. 可选分类注释       | `8-SINTAX/`（若选）  | `otus_sintax.txt`                          |
| 日志                  | —                    | `amplicon_processing.log`                  |

######################################################################################################################################################################################################################

### Amplicaon_processing2.py 脚本功能与用途
## 脚本功能及用途  
该脚本用于自动化处理扩增子（amplicon）测序数据，从原始双端 FASTQ 文件到最终 OTU 表，并可选地进行分类注释。主要流程包括：  
1. 样本解复用（demultiplex）  
2. 双端合并（merge paired‐ends）  
3. 质量过滤（quality filtering）  
4. 序列去重复（dereplication）  
5. OTU/ASV 聚类（clustering）  
6. 嵌合体检测（chimera detection）  
7. OTU 表生成（OTU table）  
8. 可选 SINTAX 分类注释（taxonomy assignment）  
此外，还提供了：  
- OTU 序列重新标记与可选二次聚类  
- 稀释曲线绘制并输出 PDF  
- 低丰度 OTU 过滤并输出过滤后序列和列表

---

## 输入

### 一、命令行参数  
| 参数                | 说明                                                         | 示例                   |
|---------------------|--------------------------------------------------------------|------------------------|
| `-i, --input_dir`   | 原始双端 FASTQ 文件所在目录（绝对或相对路径）。             | `raw_data/`            |
| `-m, --metadata_file` | 元数据文件路径（TSV 格式，无表头，包含引物和文件信息）。      | `metadata.tsv`         |
| `-o, --output_dir`  | 结果输出目录（若不存在则自动创建）。                          | `results/`             |
| `-t, --threads`     | 并行线程数（整数）。                                          | `4`                    |

### 二、元数据文件格式（TSV，无表头）  
脚本假设元数据按以下列顺序排列：  
- **run_id**：测序运行编号  
- **sample_id**：样本编号，用于命名中间及输出文件  
- **forward_primer**：正向引物序列  
- **reverse_primer**：反向引物序列  
- **forward_file**：正向 FASTQ 文件名（相对于 `input_dir`）  
- **reverse_file**：反向 FASTQ 文件名（相对于 `input_dir`）

#### 示例内容  
```tsv
run01    SampleA    AGAGTTTGATCMTGGCTCAG    CGGTTACCTTGTTACGACTT    SampleA_R1.fastq.gz    SampleA_R2.fastq.gz  
run01    SampleB    AGAGTTTGATCMTGGCTCAG    CGGTTACCTTGTTACGACTT    SampleB_R1.fastq.gz    SampleB_R2.fastq.gz  
``` 

### 三、交互式参数  
脚本在运行过程中还会交互式提示如下参数：  
- **PCR 产物长度类型**（范围或固定值）及相应长度数值  
- **聚类策略**：1. UPARSE3（OTU）  2. UNOISE3（ASV）  
- **嵌合体检测方法**：1. de novo  2. 参考数据库（需提供数据库路径）  
- **是否执行 SINTAX 分类注释**（yes/no；若 yes 则输入分类数据库路径）  
- **是否进行第二次 UPARSE3 聚类**（若选择，则输入相似性阈值）  
- **稀释曲线绘制**（自动读取样本 ID 并调用外部 R 脚本）  
- **低丰度 OTU 过滤**：最低计数阈值、最低频率阈值

---

## 输出

脚本在指定的 `output_dir` 下，按照处理步骤生成子目录和文件：

| 步骤               | 子目录            | 主要输出文件                                 |
|--------------------|-------------------|----------------------------------------------|
| 1. 解复用          | `1-demultiplex/`  | `SampleA.R1.fastq`, `SampleA.R2.fastq`       |
| 2. 合并            | `2-merge/`        | `SampleA.merged.fastq`                       |
| 3. 质量过滤        | `3-quality/`      | `SampleA.filtered.fasta`                     |
| 4. 去重复          | `4-dereplicate/`  | `SampleA.derep.fasta`, `all_samples_derep.fasta` |
| 5. 聚类            | `5-cluster/`      | `otus.fasta` (或 `centroids.fasta`)          |
| 6. 嵌合体检测      | `6-chimera/`      | `otus_nochim.fasta`                          |
| 7. OTU 表生成      | `7-OTU/`          | `otu_table.txt`                              |
| 8. 可选分类注释    | `8-SINTAX/`       | `otus_sintax.txt`                            |
| 新功能–二次聚类    | `7-OTU/`          | `otus.fasta`, `otutab.txt`                   |
| 新功能–稀释曲线    | `7-OTU/`          | `rarefaction_curve.pdf`, `rare.txt`          |
| 新功能–低丰度过滤  | `7-OTU/`          | `otutab.filter.txt`, `otus.filter.fasta`, `list.filter` |
| 日志               | —                 | `amplicon_processing.log`                    |

## 脚本依赖

- Python 3  
- 外部工具：`cutadapt`, `vsearch`, `usearch`, `seqkit`, `csvtk`  
- Python 库：`pandas`, `argparse`, `logging`
######################################################################################################################################################################################################################
### Sample Data
barcoding_corrected.txt    元数据     
id.sample    样本ID(元数据的第二列)     
1.fq.gz     正测序数据     
2.fq.gz     反测序数据     
sintax.zip    注释数据库
