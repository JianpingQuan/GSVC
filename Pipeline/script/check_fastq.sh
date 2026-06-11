#!/bin/bash
# 递归检查 fastq.gz 文件完整性，并收集坏掉的样本 ID 到 list.txt

folder=$1

if [ ! -d "$folder" ]; then
    echo "❌ 输入路径不是目录: $folder"
    exit 1
fi

bad_list="list.txt"
> "$bad_list"   # 清空旧的 list.txt

# 遍历所有 fastq.gz 文件
find "$folder" -type f -name "*.fastq.gz" | while read -r f; do
    echo "Checking $f ..."
    sample=$(basename "$f" .fastq.gz)  # 去掉后缀

    if gzip -t "$f" 2>/dev/null; then
        lines=$(zcat "$f" | wc -l)
        if (( lines % 4 == 0 )); then
            echo "[OK]   $f  ($lines lines)"
        else
            echo "[WARN] $f  ($lines lines, not multiple of 4)"
            echo "$sample" >> "$bad_list"
        fi
    else
        echo "[BAD]  $f  (gzip error)"
        echo "$sample" >> "$bad_list"
    fi
done

echo "✅ 检测完成，坏掉的样本 ID 已保存到 $bad_list"
