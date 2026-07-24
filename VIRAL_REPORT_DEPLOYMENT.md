# 病毒多样性报告模块：Linux 部署与运行

本模块建立在已有的 `contig-ui` Conda 环境和本项目 `.venv` 网页环境之上。网页端仍由 `.venv/bin/streamlit` 启动；病毒分析命令由启动 Streamlit 时所在的 Conda 环境提供。因此，推荐将生物信息学工具安装到**启动网页服务的同一个 Conda 环境**中。

## 1. 更新软件目录

将本次升级后的完整 `contig_pipeline` 目录同步到 Linux 服务器。同步后进入软件目录并赋予新增脚本可执行权限：

```bash
cd /home/zhaowl/Software/contig_pipeline
chmod u+x install.sh scripts/*.sh
PYTHON_BIN=python bash install.sh
```

不要删除已有的 assembly、cleandata 或旧 `.contig_pipeline` 运行记录。

## 2. 安装病毒分析工具

以 fat2 为例，先进入原先启动网页服务的环境：

```bash
conda activate contig-ui
cd /home/zhaowl/Software/contig_pipeline
```

如果服务器已安装 `mamba`，推荐：

```bash
mamba env update -n contig-ui -f deploy/contig-viral-tools.yml
```

没有 `mamba` 时可改用：

```bash
conda env update -n contig-ui -f deploy/contig-viral-tools.yml
```

验证：

```bash
for tool in genomad checkv coverm vclust minimap2 samtools; do
  command -v "$tool"
  "$tool" --help >/dev/null 2>&1 || echo "请检查：$tool"
done
```

默认报告不要求 VirSorter2。只有在 UI 中勾选“VirSorter2 交叉验证”时才安装它，并按其官方说明完成数据库部署：

```bash
mamba install -n contig-ui -c conda-forge -c bioconda virsorter=2
virsorter setup -d /home/zhaowl/Database/virsorter2 -j 8
```

如果 VirSorter2 2.2.4 在新版 Mamba 中报出 `Non-conda folder exists at prefix`，不要反复在同一数据库目录执行 setup；改用独立 VirSorter2 环境与 `--skip-deps-install` 方案，见 [VIRSORTER2_ISOLATED_INSTALL.md](VIRSORTER2_ISOLATED_INSTALL.md)。

## 3. 下载并配置数据库

数据库只需管理员下载一次；建议放在所有运行账号可读、普通任务账号不可写的位置。

```bash
mkdir -p /home/zhaowl/Database/genomad /home/zhaowl/Database/checkv
genomad download-database /home/zhaowl/Database/genomad
checkv download_database /home/zhaowl/Database/checkv
```

确认生成的实际目录名称后，编辑配置：

```bash
vi config/pipeline.env
```

除原有设置外，至少填写：

```bash
ALLOWED_DATA_ROOTS=/home/zhaowl/Work/WulabNGSData:/home/zhaowl/Database
GENOMAD_DB=/home/zhaowl/Database/genomad/genomad_db
CHECKV_DB=/home/zhaowl/Database/checkv/checkv-db-v1.5

MAX_THREADS_PER_VIRAL_TOOL=32
VIRAL_MIN_CONTIG_LEN=1000
VOTU_ANI=95
VOTU_ALIGNED_FRACTION=85
COVERM_MIN_READ_PERCENT_IDENTITY=95
COVERM_MIN_READ_ALIGNED_PERCENT=75
COVERM_MIN_COVERED_FRACTION=10
```

`GENOMAD_DB` 与 `CHECKV_DB` 必须改为下载后实际存在的目录；使用下列命令确认：

```bash
ls -ld "$GENOMAD_DB" "$CHECKV_DB"
```

若 `ALLOWED_DATA_ROOTS` 中不想加入整个数据库目录，也可以只保留数据项目根目录；数据库路径由受信任的 `pipeline.env` 配置读取，不来自网页用户输入。

## 4. 先在命令行做一次预检查

以下示例为 PE、每个样本一个 cleandata 子目录。输出路径应是一个新的空目录：

```bash
conda activate contig-ui
cd /home/zhaowl/Software/contig_pipeline

scripts/run_viral_report.sh \
  --assembly-dir /home/zhaowl/Work/WulabNGSData/project/assembly \
  --cleandata-dir /home/zhaowl/Work/WulabNGSData/project/cleandata \
  --clean-layout sample_subdirs \
  --read-type pe \
  --output-dir /home/zhaowl/Work/WulabNGSData/project/viral_report \
  --threads 8 \
  --min-contig-length 1000
```

单端数据时将 `--read-type se`；若 clean reads 平铺存放，将 `--clean-layout flat`。可选分组表是 UTF-8 TSV：

```text
sample_id	group
S01	control
S02	control
S03	treatment
```

成功后打开：

```text
<viral_report>/reports/virome_dashboard.html
```

## 5. 在网页端运行

完成上述安装后，重启 fat2 的 Streamlit 服务，使其继承 `contig-ui` 环境的 PATH：

```bash
conda activate contig-ui
cd /home/zhaowl/Software/contig_pipeline

kill "$(cat logs/streamlit-fat2-8502.pid)" 2>/dev/null || true
nohup env PATH="$CONDA_PREFIX/bin:$PATH" LANG=C LC_ALL=C \
  ./.venv/bin/streamlit run app.py \
  --server.address 127.0.0.1 \
  --server.port 8502 \
  --server.headless true \
  > logs/streamlit-fat2-8502.log 2>&1 &
echo $! > logs/streamlit-fat2-8502.pid
```

在 Windows 保持隧道：

```powershell
ssh -N -L 8502:127.0.0.1:8502 zhaowl@<fat2-IP>
```

浏览器打开 `http://127.0.0.1:8502`，选择“基于 contigs + cleandata 生成病毒多样性解读报告”。浏览器刷新不会终止任务；状态和日志保存在 `<viral_report>/.contig_pipeline/runs/`。

## 6. 资源和结果解释

- 首次运行应从 `--threads 8` 开始；geNomad、CheckV、CoverM 均会使用内存和临时磁盘。
- 不要同时对相同 `viral_report` 输出路径启动多个任务；脚本会加锁。
- `--resume` 仅跳过完整步骤；不完整的步骤输出不会自动删除或覆盖。
- 默认 vOTU 阈值为 ANI 95%、aligned fraction 85%。这是用于本系统内部一致比较的定义；报告必须同时查看 CheckV 质量、覆盖度和未分类比例。
- HTML 报告描述的是序列与相对丰度证据，不证明病毒具有感染性、活性或因果关系。
