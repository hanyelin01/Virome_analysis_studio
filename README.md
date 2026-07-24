# 高通量数据 Contig 与病毒多样性 Pipeline

该项目在 Linux 服务器运行，由 Windows 通过 SSH 隧道访问中文 Streamlit 界面。

支持：

- `rawdata → cleandata`：双端 PE fastp；支持“样本子文件夹”或“平铺 FASTQ”。
- `cleandata → contigs`：PE/SE MEGAHIT；支持“样本子文件夹”或“平铺 FASTQ”。
- `rawdata → cleandata → contigs`：PE 完整流程；fastp 的标准输出会自动交给 MEGAHIT。
- `contigs + cleandata → 全局病毒目录报告（v2）`：geNomad、VirSorter2、DIAMOND-NR-virus 三工具发现；完全去冗余 VC 目录；完整 NR/TaxonKit 分类、DAA/RMA6、全局 CheckV、ICTV 本地库精细注释、样本分发与 CoverM 回贴。详见 [VIROME_CATALOGUE_V2.md](VIROME_CATALOGUE_V2.md)。
- `DIAMOND 精细注释`：使用已有 `viral_report` 中 CheckV 筛选的候选 contigs，或用户提供的单个/多个 FASTA；支持完整 NR 库的病毒 TaxID 10239 过滤、无分类过滤或自定义 TaxID，并可生成 DAA、MEGAN RMA6、DIAMOND outfmt 6 和 TaxonKit LCA 注释表。

fastp 去接头默认使用双端重叠与 PE 自动识别；网页也可在实验记录明确时选择版本化的建库接头方案。每个样本会保存 fastp 报告、识别序列、目录匹配、来源等级和接头目录哈希。接头目录的审核、维护和回退规则见 [ADAPTER_CATALOG_MAINTENANCE.md](ADAPTER_CATALOG_MAINTENANCE.md)。

## 安装

将该目录复制到 Linux（例如 `/opt/contig_pipeline`）。安装前必须使用 Python 3.8+（推荐 Conda Python 3.10）；完整故障处理与部署命令见 [DEPLOYMENT.md](DEPLOYMENT.md)。执行：

```bash
cd /opt/contig_pipeline
PYTHON_BIN=python bash install.sh
```

安装前后都必须编辑 `config/pipeline.env`，至少设置可操作的数据根目录：

```bash
ALLOWED_DATA_ROOTS=/srv/contig_projects:/home/hanyl/Projects
```

运行网页服务的 Linux 账号还需要在 `PATH` 中有：`fastp`、`megahit`、`bash`、`find`、`realpath` 和 `flock`。全局病毒目录模块还需要 `genomad`、`virsorter`、`checkv`、`coverm`、`diamond`、`taxonkit`，并配置 geNomad、CheckV、完整 NR、MEGAN 和本地 ICTV 参考库。

启动服务：

```bash
/opt/contig_pipeline/.venv/bin/streamlit run /opt/contig_pipeline/app.py \
  --server.address 127.0.0.1 --server.port 8501
```

Windows PowerShell：

```powershell
ssh -N -L 8501:127.0.0.1:8501 username@linux-server
```

浏览器打开 `http://127.0.0.1:8501`。

如需在服务器重启后自动运行，可编辑 `deploy/contig-pipeline.service.example` 中的服务账号与路径，复制为 `/etc/systemd/system/contig-pipeline.service`，然后执行：

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now contig-pipeline
```

服务账号应是低权限账号，例如 `task-runner`，并且它必须对配置的数据根目录具有所需的读写权限；不要以 root 运行网页服务。

## 输出和恢复规则

完整流程的 fastp 输出固定为：

```text
cleandata/<sample>/<sample>_R1.clean.fq.gz
cleandata/<sample>/<sample>_R2.clean.fq.gz
```

MEGAHIT 输出固定为：

```text
assembly/<sample>/final.contigs.fa
```

日志、参数快照、样本清单和结果汇总自动保存于：

```text
<assembly 输出路径或 cleandata 路径>/.contig_pipeline/
```

`--resume` 只跳过完整输出。部分输出绝不会自动覆盖或删除，须人工移动或清理后再运行。脚本使用 `flock` 阻止同一输出位置的并发运行。

网页提交后，pipeline 在 Linux 后台运行；刷新或重新打开浏览器不会中止它。页面会从 `.contig_pipeline/runs/` 读取当前输出位置最近一次任务的持久化状态和末尾日志；点击“刷新任务状态”即可更新显示。

全局病毒目录流程的 assembly 输入严格受本次样本 manifest 限制；`--resume` 会核对输入、数据库和结果参数契约，不会复用其他批次或其他参考库的结果。历史旧版报告目录迁移说明见 [UPGRADE_MANIFEST_SCOPE_FIX.md](UPGRADE_MANIFEST_SCOPE_FIX.md)。

## 识别规则与限制

PE 识别 `_R1_001/_R2_001`、`_R1/_R2`、`_1/_2`，以及标准 clean 命名 `_R1.clean/_R2.clean`；扩展名可为 `.fq.gz` 或 `.fastq.gz`。每个样本必须恰有一对；多条 lane 应先合并或拆分。

SE 每个样本必须恰有一个 `.fq.gz` 或 `.fastq.gz` 文件。rawdata 质控目前仅支持 PE，因为所提供的 SOP fastp 参数属于双端处理。

网页服务只调用固定总控脚本，不提供任意 Shell 命令；路径会被限制在 `ALLOWED_DATA_ROOTS`。

## 全局病毒目录报告输出（v2）

报告任务的输入是已经完成的 `assembly/<sample>/final.contigs.fa` 以及对应样本的 clean reads。结果写入用户选择的 `virome_catalogue/` 目录：

```text
virome_catalogue/
├── 01_prepared_contigs/      # 统一且可追溯的 contig 标识
├── 03_candidate_catalogue/   # 三工具证据与全局 VC 完全去冗余目录
├── 04_nr_annotation/         # 完整 NR、TaxonKit LCA、DAA/RMA6 与判定
├── 05_checkv/                # 全局质量、完整度与前噬菌体处理
├── 06_ictv_refinement/       # 本地 ICTV 参考库精细比对
├── 07_final_catalogue/       # 最终 VF 病毒片段与来源关系
├── 08_sample_results/        # VF 回分发到原始样本
├── 09_abundance/             # CoverM reads 回贴
├── reports/virome_catalogue_dashboard.html
└── .contig_pipeline/runs/     # 可恢复的任务日志和参数快照
```

HTML 报告首先呈现三工具发现交集（UpSet 风格）、NR 证据判定、分类层级、巴尔的摩分类、CheckV 与 ICTV 精细注释证据，以及最终片段回分发到样本的结果。VC/VF 不是 vOTU 或病毒物种。
