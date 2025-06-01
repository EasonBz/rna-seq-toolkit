#!/bin/bash
# RNA-seq工具包数据准备脚本

echo "🧬 RNA-seq工具包数据准备"
echo "=========================="

# 创建数据目录结构
echo "📁 创建目录结构..."
mkdir -p data/{fastq,reference,reference/star_index}
mkdir -p results
mkdir -p docs

echo "✅ 目录结构创建完成"

# 显示目录结构
echo ""
echo "📋 项目目录结构:"
tree -L 3 . 2>/dev/null || find . -type d -not -path '*/\.*' | head -20

echo ""
echo "📊 数据放置说明:"
echo "================"
echo "1. FASTQ文件 -> data/fastq/"
echo "   示例: sample1_1.fastq.gz, sample1_2.fastq.gz"
echo ""
echo "2. 参考基因组 -> data/reference/genome.fa"
echo "   下载命令: wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
echo ""
echo "3. 注释文件 -> data/reference/annotation.gtf"
echo "   下载命令: wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
echo ""
echo "4. STAR索引 -> data/reference/star_index/"
echo "   生成命令: STAR --runMode genomeGenerate --genomeDir data/reference/star_index --genomeFastaFiles data/reference/genome.fa --sjdbGTFfile data/reference/annotation.gtf --sjdbOverhang 99 --runThreadN 16"

echo ""
echo "🔧 下一步操作:"
echo "============="
echo "1. 将数据文件放置到相应目录"
echo "2. 安装依赖: pip install -r requirements.txt"
echo "3. 运行测试: python tests/test_refactored_tools.py"
echo "4. 开始分析: python examples/basic_analysis.py"

echo ""
echo "📖 详细说明请查看: data/README_data.md" 