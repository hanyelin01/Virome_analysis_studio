# CheckV前噬菌体子片段ID恢复修复

## 现象与原因

CheckV对包含前噬菌体区域的输入contig会在`proviruses.fna`中输出切出的子片段。例如：

```text
父contig：PR10_HYY_BJX-03_lane1__k79_43359
子片段：  PR10_HYY_BJX-03_lane1__k79_43359_1 291-1036/1036
```

`_1`表示CheckV切出的第1个前噬菌体片段。旧版只用子片段ID精确查询`contig_provenance.tsv`，因此把它误报为“missing provenance or sample not in manifest”。

## 修复行为

- 始终优先精确匹配完整ID，避免误截断本身以`_1`结尾的原始contig。
- 精确匹配失败时，只尝试移除一个末尾`_<正整数>`。
- 只有父ID存在于provenance，而且CheckV质量表明确标记父contig为`provirus=Yes`时，才接受该恢复。
- `candidate_metadata.tsv`新增`provenance_sequence_id`和`provenance_match`，记录父ID及`checkv_provirus_suffix`恢复方式。
- 只有样本拆分完全成功且未分配数为0时才写出`split_complete.json`。失败遗留的`split_summary.tsv`不再被当作完成标志。

## 当前8样本任务

真实数据复核结果：

```text
Restored 88 CheckV candidates to 8 sample directories
unassigned_candidate_count=0
```

失败的746 bp前噬菌体子片段已正确恢复到`PR10_HYY_BJX-03_lane1`，随后按1,000 bp阈值归入短序列排除统计。

部署补丁后，可在网页保持相同输入和输出目录：

```text
/home/hanyl/Projects/0NIPBCRlabNGSData/1LBinfodata/viral_report_8samples_v2
```

勾选`resume`重新提交第四步。prepare、geNomad和CheckV会安全跳过；由于旧失败结果没有`split_complete.json`，vOTU阶段会自动重建样本拆分，再继续聚类、CoverM和报告生成。
