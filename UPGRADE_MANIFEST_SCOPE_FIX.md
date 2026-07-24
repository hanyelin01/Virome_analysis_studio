# 病毒报告样本范围与安全续跑修复

## 修复原因

旧版 `04_prepare_viral_contigs.sh` 在 assembly 根目录下扫描全部 `*/final.contigs.fa`，没有使用第四步预检生成的本次样本 manifest。若 assembly 中保留旧批次目录，而报告输出又启用 `--resume`，旧的 `01_prepared_contigs`、geNomad 和 CheckV 结果可能包含本次 manifest 以外的序列。流程会在恢复样本归属时停止并生成 `unassigned_candidates.tsv`。

## 新版行为

- 只处理当前 `sample_manifest.tsv` 明确列出的 `sample_id` 和 `assembly_dir`。
- 检查 manifest 路径必须与所选 assembly 根目录及样本名一致。
- 每次准备 contig 后保存 `01_prepared_contigs/preparation_inputs.json`。
- 指纹包含样本集合、每个 `final.contigs.fa` 的路径/大小/修改时间和最短长度参数。
- `--resume` 只有在指纹完全一致且四个准备产物完整时才允许跳过。
- 旧版没有指纹的报告目录会被安全拒绝，不再静默复用。

## 部署

```bash
rsync -av \
  --exclude '.git/' \
  --exclude '.venv/' \
  --exclude 'config/pipeline.env' \
  --exclude 'logs/' \
  --exclude '__pycache__/' \
  /home/hanyl/Work/Software/contig_pipeline.update/ \
  /home/hanyl/Work/Software/contig_pipeline/

cd /home/hanyl/Work/Software/contig_pipeline
chmod 755 scripts/run_pipeline.sh scripts/run_viral_report.sh scripts/run_fine_annotation.sh
python3 -m unittest discover -s tests -v
```

Git已记录三个网页入口脚本为可执行文件（`100755`）；部署后的 `chmod` 是对目标文件系统的额外保险。

## 当前失败任务的处理

不要删除单条未分配候选，也不要继续复用原报告目录：

```text
/home/hanyl/Projects/0NIPBCRlabNGSData/1LBinfodata/viral_report
```

保留它作为故障证据，在网页第四步选择一个新的输出目录，例如：

```text
/home/hanyl/Projects/0NIPBCRlabNGSData/1LBinfodata/viral_report_8samples_v2
```

重新提交相同8个样本。新版会从manifest重新建立`01_prepared_contigs`，assembly根目录中的其他历史样本不会进入新报告。

若必须沿用原目录，需要先完整归档旧目录，再整体移走`01_prepared_contigs`到`04_sample_votu`及报告产物；不建议在生产数据上做局部删除，因为各阶段结果具有链式依赖。
