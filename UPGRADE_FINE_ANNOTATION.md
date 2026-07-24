# DIAMOND 精细注释升级部署（Ubuntu 24.04 / nipb-x1）

本次升级增加独立的“DIAMOND 精细注释”任务。它不会重跑、删除或覆盖既有的 `cleandata`、`assembly`、`viral_report` 和历史日志。

本文假定软件目录为：

```text
/home/hanyl/Work/Software/contig_pipeline
```

## 1. 同步新版软件并保留服务器配置

先备份服务器专属配置，再同步新版完整目录。不要用示例文件覆盖已有 `config/pipeline.env`，其中保存的是服务器数据库路径。

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline

cp -a config/pipeline.env "config/pipeline.env.backup.$(date +%Y%m%d_%H%M%S)"
chmod u+x install.sh scripts/*.sh
PYTHON_BIN=python bash install.sh
```

若通过 `rsync` 同步，建议排除本机配置和网页虚拟环境：

```bash
rsync -av \
  --exclude '.venv/' \
  --exclude 'config/pipeline.env' \
  /新版/contig_pipeline/ \
  /home/hanyl/Work/Software/contig_pipeline/
```

同步完成后，重新执行上面的 `chmod` 与 `install.sh`。只有安装器提示 Python 主版本不同，才删除该软件目录内的 `.venv` 后重建；不要删除项目数据目录。

## 2. 安装新增工具

新版部署环境加入了 DIAMOND 和 TaxonKit：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline

mamba env update -n contig-ui -f deploy/contig-viral-tools.yml
# 没有 mamba 时：conda env update -n contig-ui -f deploy/contig-viral-tools.yml

command -v diamond taxonkit
diamond version
taxonkit version
```

MEGAN 的 `daa2rma` 不在此 Conda 环境中；只有需要 RMA6 文件时，才按你已有的 MEGAN 安装提供其程序和 mapping database。

## 3. 配置 NR、TaxonKit 和 MEGAN

编辑配置：

```bash
cd /home/hanyl/Work/Software/contig_pipeline
vi config/pipeline.env
```

在文件末尾添加下列条目，并将路径改成真实路径：

```bash
# 包含项目、viral_report 和自定义候选 contigs 的共同根目录。
ALLOWED_DATA_ROOTS=/home/hanyl/Work:/home/hanyl/Database

# 这是完整 NR DIAMOND 数据库，不是 virus.dmnd。
DIAMOND_NR_DB=/home/hanyl/Database/diamond/nr.dmnd

# “仅病毒”选项将在命令中加入 --taxonlist 10239。
DIAMOND_DEFAULT_TAXONLIST=10239
DIAMOND_EVALUE=1e-5
DIAMOND_NR_MAX_TARGET_SEQS=25
DIAMOND_SENSITIVITY=more-sensitive

# 文件夹中必须有 names.dmp、nodes.dmp、merged.dmp、delnodes.dmp。
TAXONKIT_DB=/home/hanyl/Database/taxonkit

# 仅生成 RMA6 时需要：
MEGAN_DAA2RMA=/home/hanyl/Software/MEGAN/tools/daa2rma
MEGAN_MAP_DB=/home/hanyl/Database/megan/megan-nr-r2.db
```

NR `.dmnd` 必须在建库时包含 NCBI taxonomy mapping；否则 DIAMOND 的 `--taxonlist` 与 `staxids` 等分类字段不可用。验证：

```bash
source config/pipeline.env
ls -lh "$DIAMOND_NR_DB"
diamond dbinfo --db "$DIAMOND_NR_DB"

ls "$TAXONKIT_DB"/names.dmp "$TAXONKIT_DB"/nodes.dmp \
   "$TAXONKIT_DB"/merged.dmp "$TAXONKIT_DB"/delnodes.dmp

# 仅 RMA6 模式需要
ls -l "$MEGAN_DAA2RMA" "$MEGAN_MAP_DB"
```

若还没有 TaxonKit taxonomy dump：

```bash
mkdir -p /home/hanyl/Database/taxonkit /tmp/taxonkit-download
cd /tmp/taxonkit-download
wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
cp names.dmp nodes.dmp merged.dmp delnodes.dmp /home/hanyl/Database/taxonkit/
```

## 4. 首次命令行测试

推荐先对已有 `viral_report` 的 CheckV 候选 contigs 测试。默认在完整 NR 中限制病毒分类号 10239：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline

scripts/run_fine_annotation.sh \
  --source checkv \
  --viral-report-dir /home/hanyl/Work/project/viral_report \
  --mode both \
  --taxon-scope virus \
  --threads 8 \
  --max-target-seqs 25 \
  --resume
```

完成后，主报告会自动刷新并增加“DIAMOND 精细注释”区块：

```text
<viral_report>/01_diamond_megan/viral_candidates.nr.daa
<viral_report>/01_diamond_megan/viral_candidates.nr.rma6
<viral_report>/02_diamond_nr_taxonomy/nr_virus_hits.outfmt6.tsv
<viral_report>/02_diamond_nr_taxonomy/contig_taxonomy_lca.tsv
<viral_report>/reports/virome_dashboard.html
```

若需要检索完整 NR，不限制分类范围：

```bash
scripts/run_fine_annotation.sh \
  --source checkv \
  --viral-report-dir /home/hanyl/Work/project/viral_report \
  --mode taxonomy \
  --taxon-scope none \
  --threads 8 \
  --max-target-seqs 25
```

自定义候选 contigs 产生独立补充报告，不会与原 vOTU 丰度矩阵混合：

```bash
scripts/run_fine_annotation.sh \
  --source custom \
  --custom-input /home/hanyl/Work/project/custom_candidates.fna \
  --custom-input-type file \
  --output-dir /home/hanyl/Work/project/custom_annotation_001 \
  --mode taxonomy \
  --taxon-scope custom \
  --taxonlist 10239,2157 \
  --threads 8 \
  --max-target-seqs 25
```

文件夹输入使用 `--custom-input-type directory`；第一层的 `.fa`、`.fna`、`.fasta` 会按文件名排序合并。重复 contig ID 会终止任务，防止产生不可靠的合并结果。

## 5. 重启网页服务

代码升级后，必须重启 Streamlit 才能显示左侧的“⑤ DIAMOND 精细注释”：

```bash
conda activate contig-ui
cd /home/hanyl/Work/Software/contig_pipeline
mkdir -p logs

if [[ -f logs/streamlit-8501.pid ]] && kill -0 "$(cat logs/streamlit-8501.pid)" 2>/dev/null; then
  kill "$(cat logs/streamlit-8501.pid)"
fi

nohup env PATH="$CONDA_PREFIX/bin:$PATH" LANG=C LC_ALL=C \
  ./.venv/bin/streamlit run app.py \
  --server.address 127.0.0.1 \
  --server.port 8501 \
  --server.headless true \
  > logs/streamlit-8501.log 2>&1 &
echo $! > logs/streamlit-8501.pid
```

Windows PowerShell：

```powershell
ssh -N -L 8501:127.0.0.1:8501 hanyl@<Linux服务器IP>
```

浏览器打开 `http://127.0.0.1:8501`，在左侧选择“⑤ DIAMOND 精细注释”。

## 6. 使用规则

- `仅病毒`：完整 NR + `--taxonlist 10239`。
- `不限制分类范围`：完整 NR，不传递 `--taxonlist`。
- `自定义 NCBI TaxID`：一个或多个正整数，以英文逗号分隔，例如 `10239,2157`。
- TaxonKit LCA 需要多个命中，默认保留 25 个；需要 LCA 时不要设置为 1。
- 同一输出目录一次只能运行一个任务。任务状态和日志保存在 `<输出目录>/.contig_pipeline/runs/`；刷新浏览器不会终止后台任务。
