# 逐样本病毒报告与批次检出总览：升级说明

本升级将旧的“全局 vOTU 丰度/多样性报告”改为实验室日常筛查更适合的模式：

```text
每个样本独立：CheckV 候选序列 → 本地 vOTU → 本样本 reads 回贴 → 样本报告
整批样本：汇总病毒科/属/种出现在哪些样本 → 可点击进入样本报告
```

批次总览不是跨样本生态统计，也不会把不同样本中的候选序列自动认定为同一病毒株。

## 1. 同步升级后的软件目录

在 Linux 服务器上，软件目录为：

```bash
cd /home/hanyl/Work/Software/contig_pipeline
```

请完整同步本次升级后的 `contig_pipeline` 目录，尤其不要遗漏以下新增文件：

```text
scripts/helpers/split_viral_candidates_by_sample.py
scripts/helpers/build_local_votu_catalogue.py
scripts/helpers/build_sample_votu_summary.py
scripts/helpers/build_viral_report.py
scripts/07_votu_abundance.sh
scripts/08_generate_viral_report.sh
scripts/run_viral_report.sh
app.py
deploy/contig-viral-tools.yml
```

同步后执行：

```bash
chmod u+x install.sh scripts/*.sh
PYTHON_BIN=python bash install.sh
```

此操作不会删除已有的 assembly、cleandata、geNomad、CheckV 或旧报告结果。

## 2. 安装 Vclust

旧版 `skani` 不再用于 vOTU 主流程。新流程使用 Vclust 完成每个样本内部的病毒候选去冗余。

```bash
conda activate contig-ui
mamba install -y -n contig-ui -c conda-forge -c bioconda vclust

vclust --help
coverm --version
```

若服务器没有 Mamba：

```bash
conda install -y -n contig-ui -c conda-forge -c bioconda vclust
```

## 3. 更新配置

将 `config/pipeline.env.example` 中新增的三项补充到现有的
`config/pipeline.env`：

```bash
VOTU_POST_CHECKV_MIN_LEN=1000
VOTU_IMPORTANCE_RELATIVE_ABUNDANCE=5
BATCH_OVERVIEW_RANK=family
```

含义：

- `VOTU_POST_CHECKV_MIN_LEN`：CheckV 输出后再次执行的候选序列长度下限。
- `VOTU_IMPORTANCE_RELATIVE_ABUNDANCE`：单样本报告中“高优先级”的相对丰度参考阈值（百分比）。
- `BATCH_OVERVIEW_RANK`：批次总览默认层级，可为 `family`、`genus` 或 `species`；推荐 `family`。

已有的 `VOTU_ANI=95`、`VOTU_ALIGNED_FRACTION=85`、CoverM 过滤阈值应保留。

## 4. 对当前已完成 geNomad / CheckV 项目的续跑

本次升级不会复用旧的 `04_votu_abundance/`。它会新建：

```text
04_sample_votu/
```

因此，当前因 CoverM/skani 停滞留下的旧 `04_votu_abundance/` 不会阻碍新流程；请保留它作为故障记录，暂时不要删除。

在原项目输出目录中执行续跑，已完成的步骤会跳过：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline

scripts/run_viral_report.sh \
  --assembly-dir /home/hanyl/Projects/0WulabNGSData/2026BatCN_NHZY_Yunnan/03.Assembly \
  --cleandata-dir /home/hanyl/Projects/0WulabNGSData/2026BatCN_NHZY_Yunnan/02.Cleandata \
  --clean-layout sample_subdirs \
  --read-type pe \
  --output-dir /home/hanyl/Projects/0WulabNGSData/2026BatCN_NHZY_Yunnan/04.Viral_report \
  --threads 20 \
  --overview-rank family \
  --resume
```

根据实际项目路径修改三处路径。启动前应确认之前停滞的 `coverm` 和 `skani` 进程已退出。

## 5. 新的输出结构

```text
04.Viral_report/
├── 01_prepared_contigs/
├── 02_genomad/
├── 03_checkv/
├── 04_sample_votu/
│   ├── sample_status.tsv
│   ├── split_summary.tsv
│   └── <sample_id>/
│       ├── 01_candidates/
│       ├── 02_vclust/
│       ├── 03_votu/
│       ├── 04_abundance/
│       └── votu_summary.tsv
└── reports/
    ├── batch_overview.html
    ├── batch_taxa_presence.tsv
    ├── sample_report_index.tsv
    └── samples/<sample_id>.html
```

应优先打开：

```text
reports/batch_overview.html
```

该页面列出本批样本中检出的病毒科（或用户选择的属/种）及其样本。点击样本名或检出矩阵中的“●”即可打开：

```text
reports/samples/<sample_id>.html
```

如果随后运行“⑤ DIAMOND 精细注释”并选择已有 `viral_report` 作为输入，报告会自动刷新：存在有效 TaxonKit LCA 的代表序列优先显示 DIAMOND/TaxonKit 注释；未获得该注释的序列继续保留 geNomad 注释。这样不会因缺少 DIAMOND 命中而丢失候选序列。

## 6. 网页端

完成同步后需要重启 Streamlit，网页任务名称将显示为：

```text
④ 逐样本病毒报告与批次总览
```

网页端可选择批次总览默认显示“病毒科、病毒属或病毒种”。推荐保持“病毒科”，并将更细层级结果作为注释线索而非绝对鉴定。
