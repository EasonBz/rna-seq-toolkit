"""
生物信息学分析通用工具模块
包含重复使用的功能，如工具检查、文件查找、MultiQC汇总等
"""

import os
import subprocess
import shutil
from typing import List, Dict, Optional, Tuple
from glob import glob
import json
from datetime import datetime


class BioToolsValidator:
    """生物信息学工具验证器"""
    
    @staticmethod
    def check_required_tools(tools: List[str]) -> List[str]:
        """
        检查必需工具是否可用
        
        参数:
            tools: 工具名称列表
            
        返回:
            List[str]: 缺失的工具列表
        """
        missing_tools = []
        for tool in tools:
            if shutil.which(tool) is None:
                missing_tools.append(tool)
        return missing_tools
    
    @staticmethod
    def validate_tools_or_raise(tools: List[str]) -> None:
        """
        验证工具可用性，如果有缺失则抛出异常
        
        参数:
            tools: 工具名称列表
            
        抛出:
            EnvironmentError: 当有工具缺失时
        """
        missing_tools = BioToolsValidator.check_required_tools(tools)
        if missing_tools:
            raise EnvironmentError(f"缺少必需工具: {', '.join(missing_tools)}。请先安装这些工具。")
    
    @staticmethod
    def check_files_exist(files: List[str]) -> List[str]:
        """
        检查文件是否存在
        
        参数:
            files: 文件路径列表
            
        返回:
            List[str]: 不存在的文件列表
        """
        missing_files = []
        for file_path in files:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
        return missing_files
    
    @staticmethod
    def validate_files_or_raise(files: List[str]) -> None:
        """
        验证文件存在性，如果有缺失则抛出异常
        
        参数:
            files: 文件路径列表
            
        抛出:
            FileNotFoundError: 当有文件缺失时
        """
        missing_files = BioToolsValidator.check_files_exist(files)
        if missing_files:
            raise FileNotFoundError(f"以下文件不存在: {', '.join(missing_files)}")


class FileFinder:
    """文件查找工具"""
    
    @staticmethod
    def find_fastq_files(input_dir: str, patterns: Optional[List[str]] = None) -> List[str]:
        """
        查找FASTQ文件
        
        参数:
            input_dir: 输入目录
            patterns: 文件模式列表，默认为常见FASTQ格式
            
        返回:
            List[str]: 找到的FASTQ文件路径列表
        """
        if patterns is None:
            patterns = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
        
        fastq_files = []
        for pattern in patterns:
            # 搜索当前目录
            fastq_files.extend(glob(os.path.join(input_dir, pattern)))
            # 搜索子目录
            fastq_files.extend(glob(os.path.join(input_dir, "*", pattern)))
        
        return sorted(fastq_files)
    
    @staticmethod
    def find_paired_fastq_files(input_dir: str) -> List[Tuple[str, str, str]]:
        """
        查找配对的FASTQ文件
        
        参数:
            input_dir: 输入目录
            
        返回:
            List[Tuple[str, str, str]]: (样本名, R1文件路径, R2文件路径) 的列表
        """
        # 查找R1文件
        r1_patterns = ["*_1.fastq.gz", "*_1.fastq", "*_R1.fastq.gz", "*_R1.fastq"]
        r1_files = []
        for pattern in r1_patterns:
            r1_files.extend(glob(os.path.join(input_dir, pattern)))
            r1_files.extend(glob(os.path.join(input_dir, "*", pattern)))
        
        paired_files = []
        for r1_file in r1_files:
            # 推断对应的R2文件
            r2_file = r1_file.replace("_1.fastq", "_2.fastq").replace("_R1.fastq", "_R2.fastq")
            
            if os.path.exists(r2_file):
                sample_name = os.path.basename(r1_file).split("_")[0]
                paired_files.append((sample_name, r1_file, r2_file))
        
        return paired_files
    
    @staticmethod
    def find_trimmed_paired_files(input_dir: str) -> List[Tuple[str, str, str]]:
        """
        查找修剪后的配对FASTQ文件
        
        参数:
            input_dir: 输入目录
            
        返回:
            List[Tuple[str, str, str]]: (样本名, R1文件路径, R2文件路径) 的列表
        """
        r1_files = glob(os.path.join(input_dir, "*_1.trimmed.fastq.gz"))
        if not r1_files:
            r1_files = glob(os.path.join(input_dir, "*_R1.trimmed.fastq.gz"))
        
        paired_files = []
        for r1_file in r1_files:
            r2_file = r1_file.replace("_1.trimmed", "_2.trimmed").replace("_R1.trimmed", "_R2.trimmed")
            
            if os.path.exists(r2_file):
                sample_name = os.path.basename(r1_file).split("_")[0]
                paired_files.append((sample_name, r1_file, r2_file))
        
        return paired_files
    
    @staticmethod
    def find_bam_files(input_dir: str, pattern: str = "*.sorted.bam") -> List[str]:
        """
        查找BAM文件
        
        参数:
            input_dir: 输入目录
            pattern: 文件模式
            
        返回:
            List[str]: 找到的BAM文件路径列表
        """
        return sorted(glob(os.path.join(input_dir, pattern)))


class MultiQCRunner:
    """MultiQC运行器"""
    
    @staticmethod
    def run_multiqc(input_dir: str, output_dir: str, force: bool = True) -> bool:
        """
        运行MultiQC生成汇总报告
        
        参数:
            input_dir: 输入目录（包含各种分析结果）
            output_dir: 输出目录
            force: 是否强制覆盖现有报告
            
        返回:
            bool: 是否成功运行
        """
        try:
            cmd = ["multiqc", input_dir, "-o", output_dir]
            if force:
                cmd.append("--force")
            
            subprocess.run(cmd, check=True)
            print("MultiQC汇总报告生成成功")
            return True
        except subprocess.CalledProcessError as e:
            print(f"警告: MultiQC运行失败 (退出码 {e.returncode})")
            return False
        except Exception as e:
            print(f"警告: MultiQC运行时发生错误: {e}")
            return False


class ProcessTracker:
    """处理过程跟踪器"""
    
    def __init__(self, process_name: str):
        self.process_name = process_name
        self.processed_samples = []
        self.failed_samples = []
        self.skipped_samples = []
        self.start_time = datetime.now()
    
    def add_processed(self, sample_name: str):
        """添加成功处理的样本"""
        self.processed_samples.append(sample_name)
    
    def add_failed(self, sample_name: str):
        """添加失败的样本"""
        self.failed_samples.append(sample_name)
    
    def add_skipped(self, sample_name: str):
        """添加跳过的样本"""
        self.skipped_samples.append(sample_name)
    
    def get_summary(self, output_dir: str, additional_info: Optional[Dict] = None) -> str:
        """
        生成处理摘要
        
        参数:
            output_dir: 输出目录
            additional_info: 额外信息字典
            
        返回:
            str: 摘要信息
        """
        end_time = datetime.now()
        duration = end_time - self.start_time
        
        summary = f"{self.process_name}完成\n"
        summary += f"处理时间: {duration.total_seconds():.1f}秒\n"
        summary += f"成功处理样本数: {len(self.processed_samples)}\n"
        summary += f"失败样本数: {len(self.failed_samples)}\n"
        summary += f"跳过样本数: {len(self.skipped_samples)}\n"
        summary += f"输出目录: {os.path.abspath(output_dir)}\n"
        
        if additional_info:
            for key, value in additional_info.items():
                summary += f"{key}: {value}\n"
        
        if self.failed_samples:
            summary += f"失败样本: {', '.join(self.failed_samples)}\n"
        
        if self.skipped_samples:
            summary += f"跳过样本: {', '.join(self.skipped_samples)}\n"
        
        return summary


class DirectoryManager:
    """目录管理器"""
    
    @staticmethod
    def ensure_directory(directory: str) -> str:
        """
        确保目录存在，如果不存在则创建
        
        参数:
            directory: 目录路径
            
        返回:
            str: 绝对路径
        """
        os.makedirs(directory, exist_ok=True)
        return os.path.abspath(directory)
    
    @staticmethod
    def clean_directory(directory: str) -> None:
        """
        清空目录内容
        
        参数:
            directory: 目录路径
        """
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)


class CommandRunner:
    """命令运行器"""
    
    @staticmethod
    def run_command_with_logging(
        cmd: List[str], 
        log_file: Optional[str] = None,
        capture_output: bool = False
    ) -> subprocess.CompletedProcess:
        """
        运行命令并记录日志
        
        参数:
            cmd: 命令列表
            log_file: 日志文件路径
            capture_output: 是否捕获输出
            
        返回:
            subprocess.CompletedProcess: 命令执行结果
        """
        if log_file:
            with open(log_file, "w") as log_fh:
                return subprocess.run(cmd, stderr=log_fh, check=True)
        else:
            return subprocess.run(cmd, capture_output=capture_output, check=True)


# 常用的工具组合
COMMON_TOOL_SETS = {
    "qc": ["fastqc", "multiqc"],
    "trimming": ["cutadapt", "multiqc"],
    "star_alignment": ["STAR", "samtools", "multiqc"],
    "hisat2_alignment": ["hisat2", "hisat2-build", "samtools", "multiqc"],
    "quantification_featurecounts": ["featureCounts", "multiqc"],
    "quantification_htseq": ["htseq-count", "multiqc"]
}

# 默认适配器序列
DEFAULT_ADAPTERS = {
    "truseq": [
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  # TruSeq Adapter, Index 1
        "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq Adapter, Index 2
    ],
    "nextera": [
        "CTGTCTCTTATACACATCT",  # Nextera transposase sequence
        "AGATGTGTATAAGAGACAG"   # Nextera transposase sequence (reverse)
    ]
} 