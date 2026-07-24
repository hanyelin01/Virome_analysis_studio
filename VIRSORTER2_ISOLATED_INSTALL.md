# VirSorter2 2.2.4 独立环境安装

适用于 VirSorter2 `setup` 在新版 Mamba 中报错：

```text
libmamba Non-conda folder exists at prefix
```

VirSorter2 2.2.4 会让 Snakemake 在数据库目录下创建内部 Conda 子环境；新版 Mamba 对已存在 prefix 的处理与其不兼容。此方案预装依赖到独立环境，并禁止 VirSorter2 创建内部环境。

## 1. 保留失败目录作为备份

```bash
cd /home/hanyl/Database

mv -n virsorter2 "virsorter2.failed.$(date +%Y%m%d_%H%M%S)"
mv -n virsorter2_20260723 "virsorter2_20260723.failed.$(date +%Y%m%d_%H%M%S)"
```

目录不存在时的 `mv` 报错可忽略。新安装成功后再按需清理备份。

## 2. 创建独立 VirSorter2 环境

不要将旧版本依赖混入 Streamlit 使用的 `contig-ui` 环境。

```bash
mamba create -y -n virsorter2 \
  -c conda-forge -c bioconda \
  python=3.8 virsorter=2 \
  click last ncbi-genome-download ruamel.yaml \
  prodigal=2.6 screed=1 'hmmer!=3.3.1' \
  scikit-learn=0.22.1 imbalanced-learn pandas=1.2 seaborn 'numpy<1.24'

conda activate virsorter2
virsorter --version
```

## 3. 下载数据库，但跳过内部环境创建

```bash
VIRSORTER_DB=/home/hanyl/Database/virsorter2_db
test ! -e "$VIRSORTER_DB" || { echo "Database directory already exists: $VIRSORTER_DB"; exit 1; }

virsorter setup -d "$VIRSORTER_DB" -j 8 --skip-deps-install
virsorter config --show-source
```

`--skip-deps-install` 是安全的，因为第二步已经安装了 VirSorter2 所需依赖。

## 4. 连接到 Contig Pipeline

同步本次升级后的 `contig_pipeline` 目录后，编辑：

```bash
vi /home/hanyl/Work/Software/contig_pipeline/config/pipeline.env
```

加入：

```bash
VIRSORTER_COMMAND=/home/hanyl/miniconda3/envs/virsorter2/bin/virsorter
VIRSORTER_USE_CONDA_OFF=1
```

重启 Streamlit 后，网页勾选“VirSorter2 交叉验证”即可。主流程会使用这个绝对路径并传递 `--use-conda-off`，避免再次触发内部 Mamba 环境创建。
