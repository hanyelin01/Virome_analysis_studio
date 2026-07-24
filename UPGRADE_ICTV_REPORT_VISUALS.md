# ICTV 参考字典与病毒报告视觉升级

本升级将批次页和单样本页改为“逐样本证据导航”风格；它不把不同样本的本地 vOTU 合并为同一株病毒，也不进行跨样本丰度比较。

## 新增的可审计参考层

- `config/ictv_family_genome_reference.tsv`：由 ICTV Master Species List MSL41.v1 (2025) 的 `Family` 与 `Genome` 字段生成的 427 个科级记录。
- `config/ictv_family_genome_reference.metadata.json`：记录来源文件、SHA-256、版本、生成规则和条目数。
- `config/priority_review_taxa.tsv`：可编辑的“关注科”导航清单。默认包含 `Coronaviridae`、`Paramyxoviridae`、`Orthomyxoviridae` 和 `Flaviviridae`，只用于优先打开样本报告。

参考表中的 `MIXED` / `review_required` 表示同一病毒科在 ICTV MSL 的种级记录中跨越了多个基因组组；报告显示“家族内异质”，不会从 contig 自动选择其中一个类型。

## 页面变化

### 批次总览

1. **主要病毒科**：展示全部分类单元。横条为该科对应本地 vOTU 的候选 contig 成员数，气泡为 mapped reads；旧结果没有 `count` 时会明确回退为相对丰度汇总。
2. **优先查看样本**：所有样本均保留在滚动列表中。排序依次为“关注科”检出、mapped reads、CheckV 质量支持与检出数。主行仅显示紧凑标签，例如 `[关注科 ×2]`；解释位于可折叠区域。
3. **病毒科 × 样本矩阵**：展示全部分类单元，并按 `DNA`、`RNA`、`RT`、`家族内异质`、`未归类` 分组。绿色方格仍仅表示本地 vOTU 有 reads 回贴支持。

### 单样本页

1. **Top 本地 vOTU 点线图**：位置表示相对丰度，气泡面积表示 mapped reads，颜色表示 ICTV 基因组组，金色外圈表示“关注科”。悬停可读取病毒科、属、分类来源和 CheckV 标签。
2. **完整度与覆盖证据**：横轴为 CheckV completeness（0–100%），纵轴为 `covered_bases / representative_length`（0–100%），气泡面积为 mapped reads，颜色为 CheckV 质量。两轴均为实际数值，不再使用无法直接辨认的对数覆盖轴。

## 部署更新

在 Linux 服务器完成 `rsync` 更新后，进入软件目录并执行：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline
./.venv/bin/pip install -r requirements.txt
```

`openpyxl` 仅用于将来刷新 ICTV 参考字典；日常生成报告不依赖它。现有 `config/pipeline.env` 被保留时，程序会自动使用软件目录中自带的两份参考 TSV，无需立即编辑配置。

若希望为既有结果生成新版报告：

```bash
bash scripts/08_generate_viral_report.sh \
  --output-dir /path/to/04.Viral_report \
  --manifest /path/to/04.Viral_report/.contig_pipeline/runs/<run_id>/sample_manifest.tsv \
  --refresh
```

若既有 `votu_summary.tsv` 没有 `read_count` 列，要使气泡真正表示 mapped reads，请先针对相同的报告目录重新运行步骤 07（不要使用 `--resume`），再运行上面的步骤 08。步骤 07 会重新做样本内 vOTU/回贴，不会重跑 geNomad、CheckV 或 DIAMOND。

## 更新 ICTV MSL 参考

下载 ICTV 最新的 `.xlsx` Master Species List 后，运行：

```bash
./.venv/bin/python scripts/helpers/build_ictv_family_reference.py \
  --input /path/to/ICTV_Master_Species_List.xlsx \
  --output config/ictv_family_genome_reference.tsv \
  --metadata config/ictv_family_genome_reference.metadata.json
```

随后使用步骤 08 的 `--refresh` 重新生成报告。建议把下载的 MSL 文件与输出的 metadata JSON 一同归档，以保存引用版本与校验和。

## 维护关注科清单

编辑 `config/priority_review_taxa.tsv`；只需保留以下字段：

```text
taxon_name    taxon_rank    enabled    display_badge    review_order
```

`enabled` 为 `yes` 的记录才参与排序；`display_badge` 默认建议使用简洁的“关注科”。不要把该清单解释为风险、宿主范围或致病性判断。
