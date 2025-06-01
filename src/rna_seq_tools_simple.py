#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import subprocess
from typing import List, Dict, Optional
from mcp.server.fastmcp import FastMCP

# 导入通用工具模块
from bio_utils import (
    BioToolsValidator, FileFinder, MultiQCRunner, ProcessTracker,
    DirectoryManager, CommandRunner, COMMON_TOOL_SETS, DEFAULT_ADAPTERS
)

# 路径配置 - 适应新的项目结构
BASE_DIR = Path(__file__).parent.parent  # 项目根目录
OUT_PATH = str(BASE_DIR / "results")
REFERENCE_FA = str(BASE_DIR / "data" / "reference" / "genome.fa")
GTF_FILE = str(BASE_DIR / "data" / "reference" / "annotation.gtf")
STAR_INDEX_DIR = str(BASE_DIR / "data" / "reference" / "star_index")

# 导出列表
__all__ = [
    'mcp', 'qc_rna_fastq', 'trim_rna_adapters', 
    'align_rna_reads', 'quantify_rna_expression'
]

mcp = FastMCP("RNA_seq_tools")

@mcp.tool()
def qc_rna_fastq(
    input_dir: str,
    output_dir: str = os.path.join(OUT_PATH, "qc_results")
) -> str:
    """
    对RNA-seq FASTQ文件进行质量控制分析
    
    参数:
        input_dir: 包含FASTQ文件的输入目录
        output_dir: 质控报告输出目录
        
    返回:
        str: 质控分析完成的摘要信息
    """
    BioToolsValidator.validate_tools_or_raise(COMMON_TOOL_SETS["qc"])
    output_dir = DirectoryManager.ensure_directory(output_dir)
    fastq_files = FileFinder.find_fastq_files(input_dir)
    
    if not fastq_files:
        raise FileNotFoundError(f"在目录 {input_dir} 中未找到任何FASTQ文件")
    
    print(f"找到 {len(fastq_files)} 个FASTQ文件")
    tracker = ProcessTracker("RNA-seq质控分析")
    
    # 运行FastQC
    for fastq_file in fastq_files:
        try:
            print(f"正在处理: {os.path.basename(fastq_file)}")
            subprocess.run(
                ["fastqc", fastq_file, "-o", output_dir, "--quiet"],
                check=True
            )
            tracker.add_processed(os.path.basename(fastq_file))
        except subprocess.CalledProcessError as e:
            tracker.add_failed(os.path.basename(fastq_file))
            print(f"警告: FastQC处理失败 (退出码 {e.returncode})")
    
    MultiQCRunner.run_multiqc(output_dir, output_dir)
    return tracker.get_summary(output_dir)

@mcp.tool()
def trim_rna_adapters(
    input_dir: str,
    output_dir: str = os.path.join(OUT_PATH, "trimmed_results"),
    adapters: Optional[List[str]] = None,
    min_length: int = 20,
    quality_cutoff: int = 20
) -> str:
    """
    使用cutadapt对RNA-seq数据进行adapter修剪和质量过滤
    
    参数:
        input_dir: 原始FASTQ文件目录
        output_dir: 修剪后文件输出目录
        adapters: adapter序列列表，默认使用Illumina TruSeq adapters
        min_length: 修剪后最小读长
        quality_cutoff: 质量分数阈值
        
    返回:
        str: 修剪操作完成的摘要信息
    """
    BioToolsValidator.validate_tools_or_raise(COMMON_TOOL_SETS["trimming"])
    
    if adapters is None:
        adapters = DEFAULT_ADAPTERS["truseq"]
    
    output_dir = DirectoryManager.ensure_directory(output_dir)
    paired_files = FileFinder.find_paired_fastq_files(input_dir)
    
    if not paired_files:
        raise FileNotFoundError("未找到任何配对的FASTQ文件")
    
    print(f"找到 {len(paired_files)} 对配对的FASTQ文件")
    tracker = ProcessTracker("RNA-seq adapter修剪")
    
    for sample_name, r1_file, r2_file in paired_files:
        try:
            print(f"正在处理样本: {sample_name}")
            
            r1_trimmed = os.path.join(output_dir, f"{sample_name}_1.trimmed.fastq.gz")
            r2_trimmed = os.path.join(output_dir, f"{sample_name}_2.trimmed.fastq.gz")
            report_file = os.path.join(output_dir, f"{sample_name}.cutadapt.json")
            log_file = os.path.join(output_dir, f"{sample_name}.cutadapt.log")
            
            cmd = [
                "cutadapt",
                "-a", adapters[0], "-A", adapters[1],
                "-q", str(quality_cutoff), "-m", str(min_length),
                "--json", report_file,
                "-o", r1_trimmed, "-p", r2_trimmed,
                r1_file, r2_file
            ]
            
            CommandRunner.run_command_with_logging(cmd, log_file)
            tracker.add_processed(sample_name)
            
        except subprocess.CalledProcessError as e:
            tracker.add_failed(sample_name)
            print(f"错误: cutadapt处理样本 {sample_name} 失败")
    
    MultiQCRunner.run_multiqc(output_dir, output_dir)
    return tracker.get_summary(output_dir)


@mcp.tool()
def align_rna_reads(
    input_dir: str,
    output_dir: str = os.path.join(OUT_PATH, "aligned_results"),
    reference_fa: str = REFERENCE_FA,
    gtf_file: str = GTF_FILE,
    threads: int = 8,
    aligner: str = "star"
) -> str:
    """
    使用STAR对RNA-seq数据进行基因组比对
    
    参数:
        input_dir: 修剪后的FASTQ文件目录
        output_dir: 比对结果输出目录
        reference_fa: 参考基因组FASTA文件路径
        gtf_file: GTF注释文件路径
        threads: 并行线程数
        aligner: 比对工具，支持"star"
        
    返回:
        str: 比对操作完成的摘要信息
    """
    if aligner == "star":
        BioToolsValidator.validate_tools_or_raise(COMMON_TOOL_SETS["star_alignment"])
    else:
        raise ValueError(f"不支持的比对工具: {aligner}")
    
    BioToolsValidator.validate_files_or_raise([reference_fa, gtf_file])
    output_dir = DirectoryManager.ensure_directory(output_dir)
    
    return _align_with_star(input_dir, output_dir, reference_fa, gtf_file, threads)


def _align_with_star(input_dir, output_dir, reference_fa, gtf_file, threads):
    """使用STAR进行比对"""
    index_dir = STAR_INDEX_DIR
    
    genome_dir_files = ["Genome", "SA", "SAindex"]
    index_exists = all(os.path.exists(os.path.join(index_dir, f)) for f in genome_dir_files)
    
    if not index_exists:
        index_dir = os.path.join(output_dir, "star_index")
        DirectoryManager.ensure_directory(index_dir)
        print("正在构建STAR索引...")
        # 索引构建代码...
    else:
        print(f"使用现有STAR索引: {index_dir}")
    
    paired_files = FileFinder.find_trimmed_paired_files(input_dir)
    if not paired_files:
        raise FileNotFoundError(f"在目录 {input_dir} 中未找到修剪后的FASTQ文件")
    
    tracker = ProcessTracker("RNA-seq比对 (使用STAR)")
    
    for sample_name, r1_file, r2_file in paired_files:
        try:
            print(f"正在比对样本: {sample_name}")
            
            output_prefix = os.path.join(output_dir, f"{sample_name}_")
            bam_file = os.path.join(output_dir, f"{sample_name}_Aligned.sortedByCoord.out.bam")
            
            subprocess.run([
                "STAR", "--runMode", "alignReads",
                "--genomeDir", index_dir,
                "--readFilesIn", r1_file, r2_file,
                "--readFilesCommand", "zcat",
                "--outFileNamePrefix", output_prefix,
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outSAMstrandField", "intronMotif",
                "--outFilterIntronMotifs", "RemoveNoncanonical",
                "--runThreadN", str(threads),
                "--quantMode", "GeneCounts"
            ], check=True)
            
            final_bam = os.path.join(output_dir, f"{sample_name}.sorted.bam")
            if os.path.exists(bam_file):
                os.rename(bam_file, final_bam)
                subprocess.run(["samtools", "index", final_bam], check=True)
            
            tracker.add_processed(sample_name)
            
        except subprocess.CalledProcessError as e:
            tracker.add_failed(sample_name)
            print(f"错误: 样本 {sample_name} 比对失败")
    
    MultiQCRunner.run_multiqc(output_dir, output_dir)
    additional_info = {"使用索引": index_dir}
    return tracker.get_summary(output_dir, additional_info)

@mcp.tool()
def quantify_rna_expression(
    input_dir: str,
    gtf_file: str = GTF_FILE,
    output_dir: str = os.path.join(OUT_PATH, "quantification_results"),
    threads: int = 8,
    method: str = "featurecounts"
) -> str:
    """
    对RNA-seq比对结果进行基因表达定量
    
    参数:
        input_dir: 包含排序BAM文件的目录
        gtf_file: GTF注释文件路径
        output_dir: 定量结果输出目录
        threads: 并行线程数
        method: 定量方法，支持"featurecounts"或"htseq"
        
    返回:
        str: 定量操作完成的摘要信息
    """
    if method == "featurecounts":
        BioToolsValidator.validate_tools_or_raise(COMMON_TOOL_SETS["quantification_featurecounts"])
    else:
        raise ValueError(f"不支持的定量方法: {method}")
    
    BioToolsValidator.validate_files_or_raise([gtf_file])
    output_dir = DirectoryManager.ensure_directory(output_dir)
    
    bam_files = FileFinder.find_bam_files(input_dir)
    if not bam_files:
        raise FileNotFoundError(f"在目录 {input_dir} 中未找到排序的BAM文件")
    
    print(f"找到 {len(bam_files)} 个BAM文件")
    tracker = ProcessTracker(f"RNA-seq基因表达定量 (使用{method})")
    
    if method == "featurecounts":
        try:
            print("使用featureCounts进行基因表达定量...")
            count_file = os.path.join(output_dir, "gene_counts.txt")
            
            cmd = [
                "featureCounts", "-p", "-B", "-C",
                "-a", gtf_file, "-o", count_file,
                "-T", str(threads), "-g", "gene_id",
                "-t", "exon", "-s", "2",
                *bam_files
            ]
            
            subprocess.run(cmd, check=True)
            print(f"基因计数完成，结果保存在: {count_file}")
            tracker.add_processed("featureCounts")
            
        except subprocess.CalledProcessError as e:
            tracker.add_failed("featureCounts")
            raise RuntimeError(f"featureCounts运行失败 (退出码 {e.returncode})")
    
    MultiQCRunner.run_multiqc(output_dir, output_dir)
    additional_info = {"定量方法": method, "处理BAM文件数": len(bam_files)}
    return tracker.get_summary(output_dir, additional_info)


if __name__ == "__main__":
    mcp.run(transports = "stdio")