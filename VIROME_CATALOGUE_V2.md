# 全局病毒目录流程（v2）

`virome_catalogue v2` 从已有 `cleandata` 和 `assembly/<sample>/final.contigs.fa` 开始。它取代网页中的旧“逐样本 vOTU 病毒报告”入口，但不会删除旧报告或旧脚本。

```text
cleandata + assembly
  → geNomad + VirSorter2 + DIAMOND-NR-virus
  → 全局 VC 完全去冗余候选目录
  → 完整 NR / TaxonKit LCA / DAA / RMA6
  → 全局 CheckV
  → 本地 ICTV 参考库精细注释
  → 按原始来源分发到样本 + CoverM 回贴
  → 离线 HTML 报告
```

## 输出对象和目录

```text
01_prepared_contigs/       输入 contig 与样本溯源
02_genomad/                geNomad 原始调用
02b_virsorter2/            VirSorter2 原始调用
02c_diamond_virus/         NR 病毒子集 DIAMOND 原始命中
03_candidate_catalogue/    VC 目录、调用证据与来源关系
04_nr_annotation/          完整 NR、TaxonKit LCA 和病毒证据判定
05_checkv/                 全局 CheckV 结果
06_ictv_refinement/        本地 ICTV 参考库命中
07_final_catalogue/        VF 目录
08_sample_results/         VF 回分发到原始样本
09_abundance/              各样本的 CoverM 回贴
reports/                   virome_catalogue_dashboard.html
```

- `VC_...` 是一条完全去冗余的候选核酸序列。
- `VF_...` 是 CheckV 修正后的最终病毒片段；前噬菌体可产生一个 VC 的子片段。
- VC/VF 都不是 vOTU、病毒物种或病毒株。`VC_source_mapping.tsv`、`VF_catalogue.tsv` 和 `sample_fragment_presence.tsv` 保留回溯关系。

## 证据判定

三工具只负责建立“潜在病毒序列池”。完整 NR 分类后，VC 被标记为：

| 结论 | 处理 |
|---|---|
| `confirmed_viral` | NR 存在病毒支持；进入 CheckV 和 ICTV。 |
| `putative_novel_virus` | 至少两种发现方法支持，但 NR 没有可靠命中；仍进入 CheckV。 |
| `ambiguous` | 病毒发现与全 NR 非病毒证据冲突；保留审计表。 |
| `nonviral_or_insufficient` | 不进入最终病毒目录。 |

TaxonKit LCA 和 ICTV 相似性均是序列证据，不构成临床诊断。

## ICTV 本地库维护

ICTV VMR 的 accession、下载的蛋白 FASTA 和审核后的元数据 TSV 必须一起冻结。元数据至少包含：

```text
reference_id  family  genus  species  baltimore_group
```

构建本地 DIAMOND 库：

```bash
bash scripts/build_ictv_reference_db.sh \
  --metadata /data/ictv/VMR_MSL41.reviewed.tsv \
  --protein-fasta /data/ictv/VMR_MSL41.proteins.faa \
  --output-dir /data/ictv/VMR_MSL41 \
  --version VMR_MSL41.v1.20260721
```

将命令输出的 `ICTV_REFERENCE_DMND`、`ICTV_REFERENCE_METADATA` 与 `ICTV_REFERENCE_VERSION` 写入 `config/pipeline.env`。构建器记录蛋白与元数据 SHA-256；更新时创建新目录和新版本，绝不覆盖已用于分析的参考库。

### 从官方 VMR 构建正式库

推荐固定 ICTV 发布页列出的 VMR 文件，而不是直接使用来源未知的旧 FASTA。`build_ictv_vmr_reference.py` 会把 VMR 中的 GenBank accession 分批提交到 NCBI E-utilities，获取 CDS 翻译蛋白；每条蛋白都保留 VMR 记录、来源 accession、蛋白 accession、基因组类型与巴尔的摩组。网络中断后使用 `--resume`，已缓存批次不会重复下载。

```bash
VMR_DIR=/home/hanyl/Database/ICTV/VMR_MSL41.v1.20260721

.venv/bin/python scripts/build_ictv_vmr_reference.py \
  --vmr-xlsx "$VMR_DIR/VMR_MSL41.v1.20260721.xlsx" \
  --output-dir "$VMR_DIR/staging" \
  --version VMR_MSL41.v1.20260721 \
  --resume

bash scripts/build_ictv_reference_db.sh \
  --metadata "$VMR_DIR/staging/ictv_reviewed_metadata.tsv" \
  --protein-fasta "$VMR_DIR/staging/ictv_proteins.faa" \
  --output-dir "$VMR_DIR/reference" \
  --version VMR_MSL41.v1.20260721
```

构建器输出 `ictv_unresolved_accessions.tsv` 和 `ictv_vmr_retrieval_manifest.json`。前者是没有可解析 accession、源分类冲突或 NCBI 未返回 CDS 的记录；它必须随版本保存并在报告中作为参考覆盖范围说明。`unclassified` 巴尔的摩组表示 VMR 的基因组类型无法以预设规则唯一映射，并不表示该病毒未被 ICTV 分类。

在提交 v2 前，`pipeline.env` 至少还要配置 `GENOMAD_DB`、`CHECKV_DB`、`DIAMOND_NR_DB`、`TAXONKIT_DB`、`MEGAN_DAA2RMA`、`MEGAN_MAP_DB`；运行账号的 `PATH` 必须有 `genomad`、`virsorter`、`checkv`、`diamond`、`taxonkit` 和 `coverm`。DIAMOND 必须为 2.2.4 或更高版本，以保留 `slineages` 和各分类等级字段。

### DIAMOND 性能参数

网页中可针对每次运行选择 `--threads`、`--block-size`、`--index-chunks` 与临时目录 `-t`；相同参数会传给病毒发现、完整 NR 分类、后台 DAA/RMA6 和 ICTV 精细比对，并写入每个阶段的 `parameters.env` 和 `diamond_command.sh`。默认值可在 `pipeline.env` 中维护：

```text
DIAMOND_THREADS_PER_JOB=64
DIAMOND_BLOCK_SIZE=4.0
DIAMOND_INDEX_CHUNKS=1
DIAMOND_TMPDIR=/dev/shm
```

`/dev/shm` 为内存盘，通常可减少临时 I/O；运行前诊断会检查其可写性及相对 block size 的可用空间。若内存盘空间不足，应降低 block size 或改用高速本地 SSD。DAA 后台任务与前台 NR 分类并行时，流程会将两者线程总数限制在 `MAX_TOTAL_THREADS` 内。

## 运行前诊断

每次首次部署、更新 Conda 环境或更新任一参考库后，先执行只读诊断：

```bash
conda activate contig-ui
python3 scripts/diagnose_virome_environment.py --format text
python3 scripts/diagnose_virome_environment.py --format json --write virome_catalogue_v2_readiness.json
```

也可在网页任务“④ 全局病毒发现、分类与精细注释”的“运行前：环境与参考数据库诊断”面板中执行。诊断检查：

- geNomad、VirSorter2、CheckV、DIAMOND、TaxonKit、CoverM 和 MEGAN `daa2rma` 是否可调用；
- NR、ICTV、CheckV、geNomad、TaxonKit 与 MEGAN 参考路径及 ICTV TSV 表头；
- 会影响结果的长度、线程、TaxID、e-value、CoverM 阈值和 VirSorter2 配置。

退出码 `0` 表示可以启动（可带警告），退出码 `3` 表示有阻断性问题。诊断不会下载、修改数据库或读取测序数据；它不替代真实小批次验收。
命令行诊断会检查**当前进程的 PATH**；若未激活 `contig-ui`，缺失工具的结果只说明该终端环境未加载相应 Conda 工具。systemd 部署时，应确保服务文件中的 PATH 同样包含 `contig-ui/bin`。

## 小型回归测试契约

`tests/fixtures/virome_catalogue_v2/` 是无敏感信息的最小合成契约。它不伪造或替代真实工具运行，而是固定验证 VC 去冗余、NR 证据判定、CheckV 子片段恢复、ICTV 注释、VF 回分发和报告命名边界：

```bash
python3 -m unittest tests.test_virome_environment_diagnostics tests.test_virome_catalogue_regression_contract -v
```

变更该目录的预期结果时，必须同时说明对应的科学规则或软件修复原因。建议在每次发布前另以独立、脱敏的小批次完成真实工具和参考库的端到端验收。

## 恢复规则

首次运行创建 `.contig_pipeline/virome_catalogue_contract.env`。只有输入、数据库和结果参数一致时 `--resume` 才会跳过完整步骤；改变 ICTV、NR、输入目录或阈值时请使用新的输出目录。
