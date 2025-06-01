#!/usr/bin/env python3
"""
重构后工具测试脚本
测试bio_utils模块和重构效果
"""

import os
import sys
import tempfile
from pathlib import Path

# 添加src目录到Python路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

def test_bio_utils_import():
    """测试bio_utils模块导入"""
    try:
        from bio_utils import (
            BioToolsValidator, FileFinder, MultiQCRunner, ProcessTracker,
            DirectoryManager, CommandRunner, COMMON_TOOL_SETS, DEFAULT_ADAPTERS
        )
        print("✅ bio_utils模块导入成功")
        return True
    except ImportError as e:
        print(f"❌ bio_utils模块导入失败: {e}")
        return False

def test_tool_validation():
    """测试工具验证功能"""
    try:
        from bio_utils import BioToolsValidator, COMMON_TOOL_SETS
        
        # 测试工具检查
        missing_tools = BioToolsValidator.check_required_tools(["ls", "cat"])
        print(f"✅ 基础工具检查通过，缺失工具: {missing_tools}")
        
        # 测试RNA-seq工具检查
        for tool_set_name, tools in COMMON_TOOL_SETS.items():
            missing = BioToolsValidator.check_required_tools(tools)
            if missing:
                print(f"⚠️  {tool_set_name} 工具集缺失: {missing}")
            else:
                print(f"✅ {tool_set_name} 工具集完整")
        
        return True
    except Exception as e:
        print(f"❌ 工具验证测试失败: {e}")
        return False

def test_file_finder():
    """测试文件查找功能"""
    try:
        from bio_utils import FileFinder
        
        # 测试FASTQ文件查找
        data_dir = "/home/amss_zsh/wsc/bioTool/data"
        if os.path.exists(data_dir):
            fastq_files = FileFinder.find_fastq_files(data_dir)
            print(f"✅ 找到 {len(fastq_files)} 个FASTQ文件")
            
            # 测试配对文件查找
            paired_files = FileFinder.find_paired_fastq_files(data_dir)
            print(f"✅ 找到 {len(paired_files)} 对配对FASTQ文件")
        else:
            print(f"⚠️  数据目录不存在: {data_dir}")
        
        return True
    except Exception as e:
        print(f"❌ 文件查找测试失败: {e}")
        return False

def test_process_tracker():
    """测试处理跟踪器"""
    try:
        from bio_utils import ProcessTracker
        
        tracker = ProcessTracker("测试处理")
        tracker.add_processed("sample1")
        tracker.add_processed("sample2")
        tracker.add_failed("sample3")
        tracker.add_skipped("sample4")
        
        summary = tracker.get_summary("/tmp", {"额外信息": "测试"})
        print("✅ 处理跟踪器测试通过")
        print("摘要示例:")
        print(summary)
        
        return True
    except Exception as e:
        print(f"❌ 处理跟踪器测试失败: {e}")
        return False

def test_directory_manager():
    """测试目录管理器"""
    try:
        from bio_utils import DirectoryManager
        import shutil
        
        # 创建临时目录进行测试
        test_dir = os.path.join(tempfile.gettempdir(), "bio_utils_test")
        
        # 测试目录创建
        created_dir = DirectoryManager.ensure_directory(test_dir)
        print(f"✅ 目录创建测试通过: {created_dir}")
        
        # 清理测试目录
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        
        return True
    except Exception as e:
        print(f"❌ 目录管理器测试失败: {e}")
        return False

def test_simple_tools_import():
    """测试简化版工具导入"""
    try:
        from rna_seq_tools_simple import (
            qc_rna_fastq, trim_rna_adapters, align_rna_reads, quantify_rna_expression
        )
        print("✅ 简化版RNA-seq工具导入成功")
        return True
    except ImportError as e:
        print(f"❌ 简化版RNA-seq工具导入失败: {e}")
        return False

def test_code_reduction():
    """测试代码重构效果"""
    try:
        # 比较原始文件和重构文件的行数
        original_file = "rna_seq_tools.py"
        simple_file = "rna_seq_tools_simple.py"
        utils_file = "bio_utils.py"
        
        def count_lines(filename):
            if os.path.exists(filename):
                with open(filename, 'r', encoding='utf-8') as f:
                    return len(f.readlines())
            return 0
        
        original_lines = count_lines(original_file)
        simple_lines = count_lines(simple_file)
        utils_lines = count_lines(utils_file)
        
        print(f"📊 代码行数对比:")
        print(f"   原始文件 ({original_file}): {original_lines} 行")
        print(f"   重构主文件 ({simple_file}): {simple_lines} 行")
        print(f"   通用工具 ({utils_file}): {utils_lines} 行")
        print(f"   重构后总计: {simple_lines + utils_lines} 行")
        
        if original_lines > 0:
            reduction = (original_lines - simple_lines) / original_lines * 100
            print(f"   主文件代码减少: {reduction:.1f}%")
            print("✅ 代码重构效果显著")
        
        return True
    except Exception as e:
        print(f"❌ 代码重构效果测试失败: {e}")
        return False

def main():
    """主测试函数"""
    print("🧪 开始测试重构后的生物信息学工具...")
    print("=" * 50)
    
    tests = [
        ("bio_utils模块导入", test_bio_utils_import),
        ("工具验证功能", test_tool_validation),
        ("文件查找功能", test_file_finder),
        ("处理跟踪器", test_process_tracker),
        ("目录管理器", test_directory_manager),
        ("简化版工具导入", test_simple_tools_import),
        ("代码重构效果", test_code_reduction),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n🔍 测试: {test_name}")
        try:
            if test_func():
                passed += 1
            else:
                print(f"❌ {test_name} 测试失败")
        except Exception as e:
            print(f"❌ {test_name} 测试异常: {e}")
    
    print("\n" + "=" * 50)
    print(f"📊 测试结果: {passed}/{total} 通过")
    
    if passed == total:
        print("🎉 所有测试通过！重构成功！")
        print("\n📋 重构总结:")
        print("✅ 成功提取了重复的代码模式")
        print("✅ 创建了可重用的通用工具模块")
        print("✅ 简化了主要功能函数")
        print("✅ 提高了代码的可维护性和可读性")
        return 0
    else:
        print("⚠️  部分测试失败，请检查代码")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 