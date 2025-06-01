# RNA-seq Analysis Toolkit

一个简洁、模块化的RNA-seq数据分析工具包，专为高通量转录组测序数据处理而设计。

## 🚀 特性

- **模块化设计**: 核心功能与通用工具分离，便于维护和扩展
- **简洁高效**: 重构后代码减少65.6%，保持功能完整性
- **标准流程**: 支持完整的RNA-seq分析流程
- **质量控制**: 集成FastQC和MultiQC进行质量评估
- **灵活配置**: 支持多线程处理，可自定义参数

## 📋 功能模块

### 核心分析流程
- **质量控制** (`qc_rna_fastq`): 使用FastQC进行FASTQ文件质量评估
- **接头修剪** (`trim_rna_adapters`): 使用cutadapt进行接头去除和质量过滤
- **序列比对** (`align_rna_reads`): 使用STAR进行基因组比对
- **基因定量** (`quantify_rna_expression`): 使用featureCounts进行基因表达定量

### 通用工具模块 (`bio_utils.py`)
- 工具可用性验证
- 文件查找和配对
- 进程跟踪和日志
- 目录管理
- 错误处理

## 🛠️ 依赖工具

### 必需工具
- **STAR** (>= 2.7): RNA-seq比对工具
- **samtools** (>= 1.9): BAM文件处理
- **featureCounts** (subread包): 基因表达定量
- **cutadapt** (>= 3.0): 接头修剪
- **FastQC** (>= 0.11): 质量控制
- **MultiQC** (>= 1.9): 报告汇总

### 安装依赖
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install star samtools subread fastqc

# 使用conda
conda install -c bioconda star samtools subread-featurecounts fastqc cutadapt multiqc

# 使用pip安装Python包
pip install cutadapt multiqc
```

## 📁 项目结构

```
rna-seq-toolkit/
├── src/                          # 源代码
│   ├── rna_seq_tools_simple.py   # 主要分析工具
│   └── bio_utils.py              # 通用生物信息学工具
├── tests/                        # 测试脚本
│   ├── test_rna_seq_tools.py     # 功能测试
│   └── test_refactored_tools.py  # 重构效果测试
├── data/                         # 数据目录
│   ├── fastq/                    # 原始FASTQ文件
│   ├── reference/                # 参考基因组和注释
│   │   ├── genome.fa             # 参考基因组FASTA文件
│   │   ├── annotation.gtf        # 基因注释GTF文件
│   │   └── star_index/           # STAR索引目录
│   └── README_data.md            # 数据说明文档
├── results/                      # 分析结果
├── examples/                     # 使用示例
├── docs/                         # 文档
├── requirements.txt              # Python依赖
├── .gitignore                    # Git忽略文件
└── README.md                     # 项目说明
```

## 📊 数据准备

### 1. 创建数据目录结构
```bash
mkdir -p data/{fastq,reference,reference/star_index}
```

### 2. 准备参考数据
将以下文件放置在相应目录：

```
data/reference/
├── genome.fa              # 人类基因组 (如 hg38.fa)
├── annotation.gtf         # 基因注释文件 (如 hg38.ncbiRefSeq.gtf)
└── star_index/           # STAR索引目录
    └── [STAR索引文件]
```

### 3. 准备测序数据
将FASTQ文件放置在 `data/fastq/` 目录：

```
data/fastq/
├── sample1_1.fastq.gz    # 样本1 Read1
├── sample1_2.fastq.gz    # 样本1 Read2
├── sample2_1.fastq.gz    # 样本2 Read1
├── sample2_2.fastq.gz    # 样本2 Read2
└── ...
```

## 🚀 快速开始

### 1. 克隆仓库
```bash
git clone <repository-url>
cd rna-seq-toolkit
```

### 2. 安装依赖
```bash
pip install -r requirements.txt
```

### 3. 运行测试
```bash
# 测试工具可用性和重构效果
python tests/test_refactored_tools.py

# 测试完整分析流程
python tests/test_rna_seq_tools.py
```

### 4. 运行分析
```python
import sys
sys.path.append('src')

from rna_seq_tools_simple import *

# 质量控制
qc_result = qc_rna_fastq(
    input_dir='data/fastq',
    output_dir='results/qc',
    threads=8
)

# 接头修剪
trim_result = trim_rna_adapters(
    input_dir='data/fastq',
    output_dir='results/trimmed',
    threads=8
)

# 序列比对
align_result = align_rna_reads(
    input_dir='results/trimmed',
    output_dir='results/aligned',
    threads=16
)

# 基因定量
quantify_result = quantify_rna_expression(
    input_dir='results/aligned',
    output_dir='results/quantified',
    threads=16
)
```

## 📈 性能特点

- **高比对率**: STAR比对通常可达90%+的唯一比对率
- **快速处理**: 多线程支持，充分利用计算资源
- **内存优化**: 合理的内存使用策略
- **错误处理**: 完善的错误检测和恢复机制

## 🔧 配置说明

### 默认路径配置
```python
# 在 rna_seq_tools_simple.py 中修改以下常量
OUT_PATH = "results/"                           # 输出目录
REFERENCE_FA = "data/reference/genome.fa"       # 参考基因组
GTF_FILE = "data/reference/annotation.gtf"     # 注释文件
```

### 线程数配置
建议根据服务器配置调整线程数：
- **小型服务器**: 4-8线程
- **中型服务器**: 8-16线程  
- **大型服务器**: 16-32线程

## 📝 许可证

MIT License

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！

## 📞 联系

如有问题或建议，请通过Issue联系。 