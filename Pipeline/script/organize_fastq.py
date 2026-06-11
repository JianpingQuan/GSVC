#!/usr/bin/env python3
import os
import shutil
import re

def organize_fastq_files(source_dir):
    """
    将配对的fastq文件按照样本ID整理到子文件夹中
    :param source_dir: 包含fastq文件的源目录
    """
    # 确保源目录存在
    if not os.path.exists(source_dir):
        raise FileNotFoundError(f"Source directory not found: {source_dir}")
    
    # 查找所有配对的fastq文件
    pattern = re.compile(r'^(ERR\d+)_([12])\.fastq\.gz$')
    file_pairs = {}
    
    # 扫描目录并识别配对文件
    for filename in os.listdir(source_dir):
        match = pattern.match(filename)
        if match:
            sample_id = match.group(1)
            pair_num = match.group(2)
            
            if sample_id not in file_pairs:
                file_pairs[sample_id] = {}
            
            file_pairs[sample_id][pair_num] = filename
    
    # 为每个样本创建子目录并移动文件
    for sample_id, files in file_pairs.items():
        # 创建子目录
        sample_dir = os.path.join(source_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        
        # 移动文件
        for pair_num, filename in files.items():
            src_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(sample_dir, filename)
            
            # 检查文件是否已存在
            if os.path.exists(dest_path):
                print(f"Warning: File already exists in destination, skipping: {dest_path}")
                continue
            
            # 移动文件
            shutil.move(src_path, dest_path)
            print(f"Moved: {src_path} -> {dest_path}")
    
    print("\nOrganization complete!")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Organize paired FASTQ files into sample subdirectories')
    parser.add_argument('source_dir', help='Directory containing the FASTQ files')
    
    args = parser.parse_args()
    
    try:
        organize_fastq_files(args.source_dir)
    except Exception as e:
        print(f"Error: {e}")
