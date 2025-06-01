#!/usr/bin/env python3
"""
RNA-seq工具测试脚本
测试比对和定量功能
"""

import os
import sys
import shutil
from pathlib import Path

# 添加src目录到Python路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

from rna_seq_tools_simple import align_rna_reads, quantify_rna_expression

def test_tool_availability():
    """测试必需工具的可用性"""
    print("=== 测试工具可用性 ===")
    
    # 测试STAR
    star_available = shutil.which("STAR") is not None
    print(f"STAR: {'✓' if star_available else '✗'}")
    
    # 测试samtools
    samtools_available = shutil.which("samtools") is not None
    print(f"samtools: {'✓' if samtools_available else '✗'}")
    
    # 测试featureCounts
    featurecounts_available = shutil.which("featureCounts") is not None
    print(f"featureCounts: {'✓' if featurecounts_available else '✗'}")
    
    # 测试multiqc
    multiqc_available = shutil.which("multiqc") is not None
    print(f"multiqc: {'✓' if multiqc_available else '✗'}")
    
    required_tools = [star_available, samtools_available, featurecounts_available]
    all_available = all(required_tools)
    
    print(f"\n必需工具状态: {'全部可用' if all_available else '部分缺失'}")
    return all_available

def test_data_availability():
    """测试数据文件的可用性"""
    print("\n=== 测试数据可用性 ===")
    
    data_dir = "/home/amss_zsh/wsc/bioTool/data"
    
    # 检查数据目录
    if not os.path.exists(data_dir):
        print(f"✗ 数据目录不存在: {data_dir}")
        return False
    
    # 检查FASTQ文件
    fastq_files = list(Path(data_dir).glob("*.fastq.gz"))
    print(f"FASTQ文件数量: {len(fastq_files)}")
    
    if len(fastq_files) < 2:
        print("✗ FASTQ文件数量不足")
        return False
    
    # 检查参考基因组
    ref_genome = "/home/amss_zsh/wsc/bioTool/data/hg38.fa"
    ref_exists = os.path.exists(ref_genome)
    print(f"参考基因组 (hg38.fa): {'✓' if ref_exists else '✗'}")
    
    # 检查GTF文件
    gtf_file = "/home/amss_zsh/wsc/bioTool/hg38.ncbiRefSeq.gtf"
    gtf_exists = os.path.exists(gtf_file)
    print(f"GTF注释文件: {'✓' if gtf_exists else '✗'}")
    
    return ref_exists and gtf_exists and len(fastq_files) >= 2

def test_align_function():
    """测试比对功能"""
    print("\n=== 测试比对功能 ===")
    
    # 直接使用原始数据目录
    data_dir = "/home/amss_zsh/wsc/bioTool/data"
    output_dir = "/home/amss_zsh/wsc/bioTool/results/test_aligned_results"
    
    print(f"输入目录: {data_dir}")
    print(f"输出目录: {output_dir}")
    
    try:
        # 清理之前的输出目录
        if os.path.exists(output_dir):
            import shutil
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        
        # 创建测试用的修剪文件，只使用1个样本的配对数据
        test_trimmed_dir = "/home/amss_zsh/wsc/bioTool/results/test_trimmed"
        if os.path.exists(test_trimmed_dir):
            shutil.rmtree(test_trimmed_dir)
        os.makedirs(test_trimmed_dir, exist_ok=True)
        
        # 获取原始FASTQ文件并选择较小的文件进行测试
        fastq_files = sorted(list(Path(data_dir).glob("*.fastq.gz")))
        print(f"找到 {len(fastq_files)} 个FASTQ文件")
        
        # 显示文件大小，选择较小的文件
        file_sizes = []
        for f in fastq_files:
            size_mb = f.stat().st_size / (1024*1024)
            file_sizes.append((f, size_mb))
            print(f"  {f.name}: {size_mb:.0f} MB")
        
        # 按大小排序，选择较小的文件
        file_sizes.sort(key=lambda x: x[1])
        
        # 选择最小的2个文件作为1个样本的配对数据
        if len(file_sizes) >= 2:
            # 尝试找到配对的文件
            selected_files = None
            for i in range(0, len(file_sizes)-1, 2):
                f1_name = file_sizes[i][0].name
                f2_name = file_sizes[i+1][0].name
                
                # 检查是否是配对文件
                if (f1_name.replace('_1.fastq.gz', '') == f2_name.replace('_2.fastq.gz', '') or
                    f1_name.replace('_R1.fastq.gz', '') == f2_name.replace('_R2.fastq.gz', '')):
                    selected_files = [file_sizes[i][0], file_sizes[i+1][0]]
                    break
            
            if not selected_files:
                # 如果没找到配对文件，就用前两个最小的文件
                selected_files = [file_sizes[0][0], file_sizes[1][0]]
            
            print(f"选择测试文件:")
            print(f"  R1: {selected_files[0].name} ({file_sizes[0][1]:.0f} MB)")
            print(f"  R2: {selected_files[1].name} ({file_sizes[1][1]:.0f} MB)")
            
            # 创建符合命名规范的链接（只有1个样本）
            sample_pairs = [
                (selected_files[0], "sample1_1.trimmed.fastq.gz"),
                (selected_files[1], "sample1_2.trimmed.fastq.gz")
            ]
            
            for src_file, dest_name in sample_pairs:
                dest_path = os.path.join(test_trimmed_dir, dest_name)
                os.symlink(str(src_file), dest_path)
                print(f"创建链接: {dest_name} -> {src_file.name}")
                
                # 简单的文件完整性检查
                try:
                    import subprocess
                    result = subprocess.run(['zcat', str(src_file)], 
                                          stdout=subprocess.PIPE, 
                                          stderr=subprocess.PIPE, 
                                          timeout=10)
                    if result.returncode == 0:
                        print(f"  ✓ {src_file.name} 文件完整性检查通过")
                    else:
                        print(f"  ✗ {src_file.name} 文件可能损坏")
                except subprocess.TimeoutExpired:
                    print(f"  ⚠ {src_file.name} 文件检查超时（文件较大）")
                except Exception as e:
                    print(f"  ⚠ {src_file.name} 文件检查异常: {e}")
        
        print("开始比对测试（只测试1个样本）...")
        result = align_rna_reads(
            input_dir=test_trimmed_dir,
            output_dir=output_dir,
            threads=16  # 使用16个线程加速处理
        )
        
        print("比对测试结果:")
        print(result)
        return True
        
    except Exception as e:
        print(f"比对测试失败: {e}")
        return False

def test_quantify_function():
    """测试定量功能"""
    print("\n=== 测试定量功能 ===")
    
    # 使用比对结果目录
    aligned_dir = "/home/amss_zsh/wsc/bioTool/results/test_aligned_results"
    output_dir = "/home/amss_zsh/wsc/bioTool/results/test_quantification_results"
    
    if not os.path.exists(aligned_dir):
        print(f"✗ 比对结果目录不存在: {aligned_dir}")
        print("请先运行比对测试")
        return False
    
    try:
        print("开始定量测试...")
        result = quantify_rna_expression(
            input_dir=aligned_dir,
            output_dir=output_dir,
            threads=16  # 使用16个线程加速处理
        )
        
        print("定量测试结果:")
        print(result)
        return True
        
    except Exception as e:
        print(f"定量测试失败: {e}")
        return False

def create_test_output_dir():
    """创建测试输出目录"""
    output_dir = "/home/amss_zsh/wsc/bioTool/results"
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def main():
    """主测试函数"""
    print("RNA-seq工具测试开始")
    print("=" * 50)
    
    # 创建输出目录
    create_test_output_dir()
    
    # 测试结果统计
    tests_passed = 0
    total_tests = 4
    
    # 1. 测试工具可用性
    if test_tool_availability():
        tests_passed += 1
        print("✓ 工具可用性测试通过")
    else:
        print("✗ 工具可用性测试失败")
    
    # 2. 测试数据可用性
    if test_data_availability():
        tests_passed += 1
        print("✓ 数据可用性测试通过")
    else:
        print("✗ 数据可用性测试失败")
        return
    
    # 3. 测试比对功能
    if test_align_function():
        tests_passed += 1
        print("✓ 比对功能测试通过")
    else:
        print("✗ 比对功能测试失败")
    
    # 4. 测试定量功能
    if test_quantify_function():
        tests_passed += 1
        print("✓ 定量功能测试通过")
    else:
        print("✗ 定量功能测试失败")
    
    # 总结
    print("\n" + "=" * 50)
    print(f"测试完成: {tests_passed}/{total_tests} 通过")
    
    if tests_passed == total_tests:
        print("🎉 所有测试通过！")
    else:
        print("⚠️  部分测试失败，请检查错误信息")

if __name__ == "__main__":
    main() 