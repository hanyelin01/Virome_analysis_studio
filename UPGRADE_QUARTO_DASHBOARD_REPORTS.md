# Quarto Dashboard 风格病毒报告升级

本次升级只替换报告层，不会重新执行 geNomad、CheckV、Vclust、CoverM 或 DIAMOND。已有的
`04_sample_votu/` 和 `02_diamond_nr_taxonomy/` 结果可直接重建新版报告。

## 1. 同步软件

在 Linux 服务器上同步升级文件时，继续保留本机配置和虚拟环境：

```bash
rsync -av \
  --exclude 'config/pipeline.env' \
  --exclude '.venv/' \
  --exclude 'logs/' \
  --exclude '__pycache__/' \
  /home/hanyl/Work/Software/contig_pipeline.update/ \
  /home/hanyl/Work/Software/contig_pipeline/

cd /home/hanyl/Work/Software/contig_pipeline
chmod u+x scripts/*.sh
```

本版报告只使用 Python 标准库，不增加 Plotly、Quarto、R 或系统级依赖；因此不需要重新安装
`.venv`。

## 2. 更新配置

将以下两行加入服务器上的 `config/pipeline.env`：

```bash
REPORT_THEME=quarto-scientific
REPORT_TOP_TAXA=0
REPORT_PRIORITY_FAMILIES=Coronaviridae,Paramyxoviridae,Orthomyxoviridae,Flaviviridae
```

`REPORT_TOP_TAXA=0` 表示批次条形图和检测矩阵显示全部分类单元，并在页面内提供滚动区域；
设置为正整数时才限制显示数量。`REPORT_PRIORITY_FAMILIES` 是“关注病毒科”的可编辑清单，
仅影响样本复核排序，不是感染性、人畜共患或致病风险判定。

## 2.1 可选：补算真实 mapped reads

旧版 `votu_summary.tsv` 只有相对丰度、平均覆盖度和覆盖碱基数。新版已让 CoverM 同时输出
`count`，因此可以在图中用气泡展示每个病毒科或本地 vOTU 的真实 mapped reads。若不执行此步，
新版报告仍可生成，但气泡会明确显示为“相对丰度合计（%）”，而非 reads 数。

以下命令只重跑每个样本的本地 Vclust/CoverM 丰度步骤；不会重跑 geNomad、CheckV、VirSorter2
或 DIAMOND。请在没有同一项目运行中的任务时执行：

```bash
REPORT_ROOT=/home/hanyl/Projects/PROJECT/04.Viral_report
MANIFEST=$(find "$REPORT_ROOT/.contig_pipeline/runs" -name sample_manifest.tsv -type f -printf '%T@\t%p\n' \
  | sort -n | tail -n 1 | cut -f2-)

bash /home/hanyl/Work/Software/contig_pipeline/scripts/07_votu_abundance.sh \
  --manifest "$MANIFEST" \
  --input "$REPORT_ROOT/03_checkv/viral_candidates_checkv.fna" \
  --output-dir "$REPORT_ROOT" \
  --threads 20
```

按服务器资源修改 `--threads`；该步骤会重新进行 reads 回贴，耗时取决于样本数和测序深度。

## 3. 仅重建报告

在已有病毒报告输出目录中执行。请将两个路径替换为真实路径：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline

bash scripts/08_generate_viral_report.sh \
  --output-dir /home/hanyl/Projects/PROJECT/04.Viral_report \
  --manifest /home/hanyl/Projects/PROJECT/04.Viral_report/.contig_pipeline/runs/RUN_ID/sample_manifest.tsv \
  --overview-rank family \
  --refresh
```

可用下列命令自动找到最近一次 manifest，再重建：

```bash
REPORT_ROOT=/home/hanyl/Projects/PROJECT/04.Viral_report
MANIFEST=$(find "$REPORT_ROOT/.contig_pipeline/runs" -name sample_manifest.tsv -type f -printf '%T@\t%p\n' \
  | sort -n | tail -n 1 | cut -f2-)

bash /home/hanyl/Work/Software/contig_pipeline/scripts/08_generate_viral_report.sh \
  --output-dir "$REPORT_ROOT" \
  --manifest "$MANIFEST" \
  --overview-rank family \
  --refresh
```

这一步不会重跑任何病毒识别或丰度计算步骤，通常只需读取已有 TSV 并生成 HTML/JSON。

## 4. 新增结果

```text
reports/
├── virome_dashboard.html       # 离线交互式主报告
├── batch_overview.html         # 兼容旧链接，自动跳转主报告
├── samples/<sample_id>.html    # 兼容旧链接，自动跳转对应样本页
└── data/
    ├── report_metadata.json
    ├── batch_dashboard_data.json
    ├── batch_summary.tsv
    ├── evidence_legend.tsv
    └── samples/<sample_id>.json
```

可直接在服务器文件系统中打开 `reports/virome_dashboard.html`；Windows 用户更推荐通过 Streamlit
的“④ 逐样本病毒报告与批次总览”任务页，点击“在网页内打开交互式报告”。也可以下载 HTML，
在无网络环境的浏览器中打开。

## 5. 证据边界

批次页只汇总“在某个样本中有 reads 支持的本地 vOTU”的分类线索。它不：

- 将不同样本的本地 vOTU 当作同一病毒株；
- 做跨样本病毒丰度或生态多样性比较；
- 将序列注释视为感染、宿主范围、致病性或活性证据。

每个图表均附带可展开的“如何解读与证据边界”说明；完整原始数值保留在
`votu_summary.tsv` 与 `reports/data/` 中。
