#!/usr/bin/env python3
"""
整合PhaTYP和BACPHLIP的噬菌体生活方式预测结果
使用阈值：BACPHLIP概率需要 >0.95 才认为是高置信度
"""

import pandas as pd
import numpy as np
import os

# 文件路径
phatyp_file = "/home/quanjp/software/PhaTYP/high_quality_combineFilter.phatyp.csv"
bacphlip_file = "/data2/qjp/BACPHLIP_input/all_vOTUs_lifestyle_predictions.tsv"
output_file = "integrated_lifestyle_predictions_thresh095.csv"

# 设置阈值
BACPHLIP_THRESHOLD = 0.95  # BACPHLIP高置信度阈值
PHATYP_THRESHOLD = 0.95     # PhaTYP也可以设个阈值，但它的分数普遍很高

# 读取PhaTYP结果
print("读取PhaTYP结果...")
phatyp_df = pd.read_csv(phatyp_file)
phatyp_df.columns = ['Contig', 'PhaTYP_pred', 'PhaTYP_score']
# 清理Contig名称（去除可能存在的空格）
phatyp_df['Contig'] = phatyp_df['Contig'].str.strip()
print(f"PhaTYP结果数量: {len(phatyp_df)}")

# 读取BACPHLIP结果
print("读取BACPHLIP结果...")
# 跳过第一行标题，直接读取数据
bacphlip_df = pd.read_csv(bacphlip_file, sep='\s+', skiprows=1, 
                          names=['Contig', 'Virulent_prob', 'Temperate_prob'])
# 清理Contig名称
bacphlip_df['Contig'] = bacphlip_df['Contig'].str.strip()

# 根据概率确定BACPHLIP的预测和置信度
bacphlip_df['BACPHLIP_pred'] = bacphlip_df.apply(
    lambda row: 'virulent' if row['Virulent_prob'] >= row['Temperate_prob'] else 'temperate', 
    axis=1
)
bacphlip_df['BACPHLIP_prob'] = bacphlip_df.apply(
    lambda row: max(row['Virulent_prob'], row['Temperate_prob']), 
    axis=1
)
bacphlip_df['BACPHLIP_high_conf'] = bacphlip_df['BACPHLIP_prob'] >= BACPHLIP_THRESHOLD

print(f"BACPHLIP结果数量: {len(bacphlip_df)}")
print(f"BACPHLIP高置信度结果(>{BACPHLIP_THRESHOLD}): {bacphlip_df['BACPHLIP_high_conf'].sum()}")

# 合并两个结果
print("\n合并结果...")
merged_df = pd.merge(phatyp_df, bacphlip_df, on='Contig', how='outer')

# 标准化预测值（确保都是小写）
merged_df['PhaTYP_pred'] = merged_df['PhaTYP_pred'].str.lower().fillna('missing')
merged_df['BACPHLIP_pred'] = merged_df['BACPHLIP_pred'].str.lower().fillna('missing')

# 定义整合规则函数（基于0.95阈值）
def integrate_predictions(row):
    phatyp_pred = row['PhaTYP_pred']
    phatyp_score = row.get('PhaTYP_score', 0)
    bacphlip_pred = row['BACPHLIP_pred']
    bacphlip_prob = row.get('BACPHLIP_prob', 0)
    bacphlip_high_conf = row.get('BACPHLIP_high_conf', False)
    
    # 处理缺失值
    if phatyp_pred == 'missing' and bacphlip_pred == 'missing':
        return 'no_prediction', 0, 'both_missing'
    elif phatyp_pred == 'missing':
        if bacphlip_high_conf:
            return bacphlip_pred, bacphlip_prob, 'only_bacphlip_high_conf'
        else:
            return bacphlip_pred, bacphlip_prob, 'only_bacphlip_low_conf'
    elif bacphlip_pred == 'missing':
        if phatyp_score >= PHATYP_THRESHOLD:
            return phatyp_pred, phatyp_score, 'only_phatyp_high_conf'
        else:
            return phatyp_pred, phatyp_score, 'only_phatyp_low_conf'
    
    # 两个工具都有结果
    if phatyp_pred == bacphlip_pred:
        # 预测一致
        confidence = max(phatyp_score, bacphlip_prob)
        if bacphlip_high_conf or phatyp_score >= PHATYP_THRESHOLD:
            return phatyp_pred, confidence, 'consistent_high_conf'
        else:
            return phatyp_pred, confidence, 'consistent_low_conf'
    else:
        # 预测不一致
        # 情况1：BACPHLIP高置信度（>0.95），以BACPHLIP为准
        if bacphlip_high_conf:
            return bacphlip_pred, bacphlip_prob, 'conflict_resolved_bacphlip_high'
        # 情况2：PhaTYP高置信度（>0.95）且BACPHLIP低置信度
        elif phatyp_score >= PHATYP_THRESHOLD and not bacphlip_high_conf:
            return phatyp_pred, phatyp_score, 'conflict_resolved_phatyp_high'
        # 情况3：两个都低置信度，但倾向选择概率更高的
        elif phatyp_score > bacphlip_prob:
            return phatyp_pred, phatyp_score, 'conflict_low_conf_prefer_phatyp'
        else:
            return bacphlip_pred, bacphlip_prob, 'conflict_low_conf_prefer_bacphlip'

# 应用整合规则
print("应用整合规则（BACPHLIP阈值 = 0.95）...")
results = merged_df.apply(lambda row: integrate_predictions(row), axis=1)
merged_df['Final_Pred'] = [r[0] for r in results]
merged_df['Final_Score'] = [r[1] for r in results]
merged_df['Decision'] = [r[2] for r in results]

# 重新排列列
columns_order = ['Contig', 'Final_Pred', 'Final_Score', 'Decision',
                 'PhaTYP_pred', 'PhaTYP_score', 
                 'BACPHLIP_pred', 'BACPHLIP_prob', 'BACPHLIP_high_conf']
available_cols = [col for col in columns_order if col in merged_df.columns]
merged_df = merged_df[available_cols]

# 保存结果
print(f"\n保存整合结果到 {output_file}...")
merged_df.to_csv(output_file, index=False)

# 输出详细统计信息
print("\n" + "="*60)
print("整合结果统计 (BACPHLIP阈值 = 0.95)")
print("="*60)

print("\n📊 最终预测分布：")
print(merged_df['Final_Pred'].value_counts())

print("\n📊 决策类型分布：")
decision_counts = merged_df['Decision'].value_counts()
for decision, count in decision_counts.items():
    print(f"  {decision}: {count}")

# 分析不一致的案例
print("\n🔍 预测不一致的详细分析：")
inconsistent = merged_df[merged_df['PhaTYP_pred'] != merged_df['BACPHLIP_pred']]
inconsistent = inconsistent[~inconsistent['PhaTYP_pred'].isin(['missing'])]
inconsistent = inconsistent[~inconsistent['BACPHLIP_pred'].isin(['missing'])]

if len(inconsistent) > 0:
    print(f"发现 {len(inconsistent)} 个预测不一致的序列")
    print("\n按照BACPHLIP置信度分组：")
    print(inconsistent['BACPHLIP_high_conf'].value_counts())
    
    # 显示需要特别注意的冲突（两个工具都高置信度但预测不同）
    serious_conflicts = inconsistent[
        (inconsistent['BACPHLIP_high_conf']) & 
        (inconsistent['PhaTYP_score'] >= PHATYP_THRESHOLD)
    ]
    if len(serious_conflicts) > 0:
        print(f"\n⚠️  警告：发现 {len(serious_conflicts)} 个严重冲突（两个工具都高置信度但预测不同）！")
        print(serious_conflicts[['Contig', 'PhaTYP_pred', 'PhaTYP_score', 
                                 'BACPHLIP_pred', 'BACPHLIP_prob']].to_string())
else:
    print("所有序列的预测都一致！")

# 保存严重冲突的序列供后续分析
if len(serious_conflicts) > 0:
    serious_conflicts.to_csv("serious_conflicts_need_review.csv", index=False)
    print(f"\n严重冲突已保存到 serious_conflicts_need_review.csv")

print(f"\n✅ 整合完成！结果保存在: {output_file}")
