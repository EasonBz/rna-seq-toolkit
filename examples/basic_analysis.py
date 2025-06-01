#!/usr/bin/env python3
"""
RNA-seq基本分析示例
演示如何使用RNA-seq工具包进行完整的分析流程
"""

import sys
import os
from pathlib import Path

# 添加src目录到Python路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rna_seq_tools_simple import (
    qc_rna_fastq, 
    trim_rna_adapters, 
    align_rna_reads, 
    quantify_rna_expression
)

def run_complete_analysis():
    """运行完整的RNA-seq分析流程"""
    
    print("🧬 开始RNA-seq分析流程")
    print("=" * 50)
    
    # 定义目录路径
    base_dir = Path(__file__).parent.parent
    fastq_dir = base_dir / "data" / "fastq"
    results_dir = base_dir / "results"
    
    # 检查输入数据
    if not fastq_dir.exists() or not list(fastq_dir.glob("*.fastq.gz")):
        print("❌ 错误: 未找到FASTQ文件")
        print(f"请将FASTQ文件放置在: {fastq_dir}")
        return False
    
    print(f"📁 输入目录: {fastq_dir}")
    print(f"📁 输出目录: {results_dir}")
    
    try:
        # 步骤1: 质量控制
        print("\n🔍 步骤1: 质量控制")
        qc_result = qc_rna_fastq(
            input_dir=str(fastq_dir),
            output_dir=str(results_dir / "01_qc"),
            threads=8
        )
        print("✅ 质量控制完成")
        print(qc_result)
        
        # 步骤2: 接头修剪
        print("\n✂️ 步骤2: 接头修剪")
        trim_result = trim_rna_adapters(
            input_dir=str(fastq_dir),
            output_dir=str(results_dir / "02_trimmed"),
            threads=8
        )
        print("✅ 接头修剪完成")
        print(trim_result)
        
        # 步骤3: 序列比对
        print("\n🎯 步骤3: 序列比对")
        align_result = align_rna_reads(
            input_dir=str(results_dir / "02_trimmed"),
            output_dir=str(results_dir / "03_aligned"),
            threads=16
        )
        print("✅ 序列比对完成")
        print(align_result)
        
        # 步骤4: 基因定量
        print("\n📊 步骤4: 基因定量")
        quantify_result = quantify_rna_expression(
            input_dir=str(results_dir / "03_aligned"),
            output_dir=str(results_dir / "04_quantified"),
            threads=16
        )
        print("✅ 基因定量完成")
        print(quantify_result)
        
        print("\n🎉 RNA-seq分析流程全部完成！")
        print(f"📋 结果文件位置: {results_dir}")
        
        return True
        
    except Exception as e:
        print(f"❌ 分析过程中出现错误: {e}")
        return False

def run_single_step_example():
    """单步分析示例"""
    
    print("\n📝 单步分析示例")
    print("-" * 30)
    
    # 只运行质量控制
    print("示例: 只运行质量控制")
    
    base_dir = Path(__file__).parent.parent
    fastq_dir = base_dir / "data" / "fastq"
    
    if fastq_dir.exists():
        try:
            result = qc_rna_fastq(
                input_dir=str(fastq_dir),
                output_dir=str(base_dir / "results" / "qc_only"),
                threads=4
            )
            print("✅ 单步质量控制完成")
            print(result)
        except Exception as e:
            print(f"❌ 单步分析失败: {e}")
    else:
        print("⚠️ 未找到FASTQ文件，跳过单步示例")

def main():
    """主函数"""
    
    print("RNA-seq分析工具包 - 使用示例")
    print("=" * 50)
    
    # 检查当前工作目录
    current_dir = Path.cwd()
    expected_dir = Path(__file__).parent.parent
    
    if current_dir != expected_dir:
        print(f"⚠️ 建议在项目根目录运行: {expected_dir}")
        print(f"当前目录: {current_dir}")
    
    # 选择运行模式
    print("\n请选择运行模式:")
    print("1. 完整分析流程")
    print("2. 单步分析示例")
    print("3. 退出")
    
    try:
        choice = input("\n请输入选择 (1-3): ").strip()
        
        if choice == "1":
            run_complete_analysis()
        elif choice == "2":
            run_single_step_example()
        elif choice == "3":
            print("👋 退出程序")
        else:
            print("❌ 无效选择")
            
    except KeyboardInterrupt:
        print("\n👋 用户中断，退出程序")
    except Exception as e:
        print(f"❌ 程序运行错误: {e}")

if __name__ == "__main__":
    main() 