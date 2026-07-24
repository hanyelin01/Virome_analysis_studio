# 参数传导审计与长度阈值修复

## 200 bp问题

真实运行快照证明页面设置的200 bp已传至contig准备阶段：

```text
parameters.env: MIN_CONTIG_LENGTH=200
preparation_inputs.json: "min_length": 200
```

但旧版在CheckV后另行读取`VOTU_POST_CHECKV_MIN_LEN=1000`，且网页没有显示该第二阈值。因此85,718条≥200 bp的contig进入geNomad，但188条CheckV候选中仅88条≥1,000 bp进入报告，100条被第二阈值排除。

新版默认勾选“CheckV后报告与vOTU沿用同一长度阈值”，设置200 bp会同时传给：

```text
网页 viral_min
→ run_viral_report.sh --min-contig-length 200
→ prepare_viral_contigs.py --min-length 200

网页 post_checkv_min
→ run_viral_report.sh --post-checkv-min-length 200
→ 07_votu_abundance.sh --min-length 200
→ split_viral_candidates_by_sample.py --min-length 200
```

真实数据按200 bp复核得到188条保留、0条短序列排除、0条未分配。

## 页面参数传导结论

| 工作流 | 页面参数 | 结论 |
|---|---|---|
| fastp | 并发样本、每样本线程、接头方案 | 已传至执行命令和运行快照 |
| MEGAHIT | 并发样本、每样本线程、最短contig | 已传至MEGAHIT命令和运行快照 |
| 病毒报告 | 分析线程 | 已传至geNomad、CheckV、vOTU/CoverM |
| 病毒报告 | geNomad输入最短长度 | 原本正常；现已在页面明确阶段 |
| 病毒报告 | CheckV后最短长度 | 原为隐藏配置；现已显式传导并默认与前段统一 |
| 病毒报告 | 分类层级 | 已传至报告生成 |
| 病毒报告 | VirSorter2开关 | 已传至可选步骤 |
| DIAMOND | 来源、模式、线程、TaxID范围、最大命中数 | 已传至MEGAN与TaxonKit分支 |

## 服务器配置参数

vOTU ANI、aligned fraction、CoverM三个过滤阈值和重要性阈值均实际用于工具命令。新版在页面展示其有效值、写入`parameters.env`，并写入`04_sample_votu/run_contract.json`。

该契约还记录manifest与CheckV候选FASTA哈希。resume时任一参数或输入变化都会明确拒绝旧结果，避免“页面已改参数但仍跳过旧输出”。

## 额外修复

命令行支持的`--groups-file`过去只被接收和记录，没有传入报告构建器。新版已贯通至`batch_summary.tsv`、`sample_report_index.tsv`和`batch_taxa_presence.tsv`的`group`列。

`MAX_VIRAL_PARALLEL`目前没有对应的并行执行模型，仍不应被理解为可生效的网页参数；病毒工具以单任务线程数控制资源。
