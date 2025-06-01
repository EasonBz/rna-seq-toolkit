# 数据目录说明

本目录用于存放RNA-seq分析所需的所有数据文件。请按照以下结构组织您的数据：

## 📁 目录结构

```
data/
├── fastq/                    # 原始测序数据
│   ├── sample1_1.fastq.gz   # 样本1的Read1文件
│   ├── sample1_2.fastq.gz   # 样本1的Read2文件
│   ├── sample2_1.fastq.gz   # 样本2的Read1文件
│   ├── sample2_2.fastq.gz   # 样本2的Read2文件
│   └── ...                  # 更多样本
├── reference/               # 参考基因组和注释文件
│   ├── genome.fa           # 参考基因组FASTA文件
│   ├── annotation.gtf      # 基因注释GTF文件
│   └── star_index/         # STAR比对索引
│       ├── Genome          # STAR索引文件
│       ├── SA              # STAR索引文件
│       ├── SAindex         # STAR索引文件
│       └── ...             # 其他STAR索引文件
└── README_data.md          # 本说明文件
```

## 📋 文件要求

### 1. FASTQ文件 (`fastq/`)

**命名规范**:
- 配对末端测序: `{样本名}_{1|2}.fastq.gz`
- 或者: `{样本名}_{R1|R2}.fastq.gz`
- 单端测序: `{样本名}.fastq.gz`

**示例**:
```
sample1_1.fastq.gz  # 样本1 Read1
sample1_2.fastq.gz  # 样本1 Read2
sample2_1.fastq.gz  # 样本2 Read1
sample2_2.fastq.gz  # 样本2 Read2
```

**格式要求**:
- 文件必须是gzip压缩格式 (`.fastq.gz`)
- 支持标准的FASTQ格式
- 质量编码: Phred+33 (Sanger格式)

### 2. 参考基因组 (`reference/genome.fa`)

**文件要求**:
- FASTA格式的参考基因组序列
- 推荐使用人类基因组 hg38/GRCh38
- 文件名建议: `hg38.fa` 或 `genome.fa`

**获取方式**:
```bash
# 从UCSC下载hg38基因组
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

### 3. 基因注释文件 (`reference/annotation.gtf`)

**文件要求**:
- GTF格式的基因注释文件
- 必须与参考基因组版本匹配
- 推荐使用NCBI RefSeq或GENCODE注释

**获取方式**:
```bash
# NCBI RefSeq注释 (推荐)
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# 或GENCODE注释
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz
```

### 4. STAR索引 (`reference/star_index/`)

**生成方法**:
```bash
# 创建STAR索引目录
mkdir -p data/reference/star_index

# 生成STAR索引 (需要大量内存，建议32GB+)
STAR --runMode genomeGenerate \
     --genomeDir data/reference/star_index \
     --genomeFastaFiles data/reference/genome.fa \
     --sjdbGTFfile data/reference/annotation.gtf \
     --sjdbOverhang 99 \
     --runThreadN 16
```

**注意事项**:
- STAR索引生成需要大量内存 (通常需要32GB+)
- 索引生成时间较长 (1-2小时)
- 生成的索引文件较大 (20-30GB)

## 🔧 数据准备脚本

创建 `prepare_data.sh` 脚本来自动化数据准备：

```bash
#!/bin/bash
# 数据准备脚本

# 创建目录结构
mkdir -p data/{fastq,reference,reference/star_index}

echo "数据目录结构已创建"
echo "请将您的数据文件放置在相应目录中："
echo "  - FASTQ文件 -> data/fastq/"
echo "  - 参考基因组 -> data/reference/genome.fa"
echo "  - 注释文件 -> data/reference/annotation.gtf"
echo "  - 然后运行STAR索引生成命令"
```

## ⚠️ 重要提醒

1. **文件大小**: RNA-seq数据文件通常很大，确保有足够的存储空间
2. **权限设置**: 确保程序有读取数据文件的权限
3. **路径配置**: 在 `src/rna_seq_tools_simple.py` 中更新文件路径
4. **备份数据**: 建议对原始数据进行备份

## 📊 数据质量检查

在开始分析前，建议检查数据质量：

```bash
# 检查FASTQ文件完整性
zcat data/fastq/sample1_1.fastq.gz | head -4

# 检查文件大小
ls -lh data/fastq/

# 检查参考文件
head data/reference/genome.fa
head data/reference/annotation.gtf
```

## 🆘 常见问题

**Q: FASTQ文件损坏怎么办？**
A: 使用 `zcat filename.fastq.gz | head` 检查文件完整性

**Q: 参考基因组版本不匹配？**
A: 确保基因组和注释文件来自同一版本 (如都是hg38)

**Q: STAR索引生成失败？**
A: 检查内存是否足够，通常需要32GB+内存

**Q: 文件路径错误？**
A: 检查 `src/rna_seq_tools_simple.py` 中的路径配置 