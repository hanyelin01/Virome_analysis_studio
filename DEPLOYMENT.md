# Linux 安装与部署指南

## 0. 本次报错的原因

报错环境使用的是系统 Python 3.6：日志中的路径为 `/usr/local/lib/python3.6/`。`streamlit==1.36.0` 要求 Python 3.8 或更新版本，因此 pip 只显示了仍兼容 Python 3.6 的旧版 Streamlit，最高到 1.10.0；这不是网络、pip 更新或权限问题。

不要继续在 Python 3.6 环境中尝试安装本项目，也不要将依赖降级到旧 Streamlit：当前界面按 Streamlit 1.36 开发。

以下命令均在 Linux 服务器执行。示例使用你的用户名和目录；请按实际 Conda 安装位置及数据路径调整。

## 1. 创建独立 Conda 环境（推荐）

先确认 Conda 可用：

```bash
conda --version
```

创建一个包含 Python 3.10、fastp 和 MEGAHIT 的独立环境。将生信软件和网页服务放在同一环境，能保证网页启动时也能找到这两个程序：

```bash
conda create -n contig-ui \
  -c conda-forge -c bioconda \
  python=3.10 fastp megahit -y

conda activate contig-ui
```

验证必须全部成功：

```bash
python --version
which python
fastp --version
megahit --version
```

其中 `python --version` 必须显示 `3.8` 或更高，建议为 `3.10.x`。若 `conda` 命令不存在，需要先请服务器管理员安装 Miniconda/Mambaforge，或使用服务器上已有的 Python 3.8+ Conda 环境。

## 2. 安装网页程序

进入项目目录。此处假定项目已经位于你的服务器目录：

```bash
cd /home/zhaowl/Software/contig_pipeline
conda activate contig-ui
```

你刚才的失败安装已创建了一个 Python 3.6 的 `.venv`。在已经确认 `python --version` 为 3.8+ 后，只删除这个项目内的旧虚拟环境：

```bash
rm -rf /home/zhaowl/Software/contig_pipeline/.venv
```

运行安装脚本；它会使用当前激活环境的 `python` 创建项目专属 `.venv`，并安装固定版本 `streamlit==1.36.0`。安装脚本还固定 `numpy==1.26.4`、`pandas==2.2.3` 和 `pyarrow==15.0.2`，它们均有 glibc 2.17+ 的预编译 wheel，避免老服务器尝试源码编译：

```bash
PYTHON_BIN=python bash install.sh
```

成功后执行以下验证：

```bash
./.venv/bin/python --version
./.venv/bin/streamlit --version
command -v fastp
command -v megahit
```

如果最后两条找不到程序，说明当前 shell 没有激活 `contig-ui`。重新运行 `conda activate contig-ui`；不要用系统 Python 3.6 启动 Streamlit。

## 3. 配置允许访问的数据目录

创建配置文件并编辑：

```bash
cp -n config/pipeline.env.example config/pipeline.env
vi config/pipeline.env
```

至少修改 `ALLOWED_DATA_ROOTS`。它是冒号分隔的 Linux 绝对路径白名单；网页内填写的 rawdata、cleandata 和 assembly 路径都必须位于其中。例如数据都在你的 home 目录下时：

```bash
ALLOWED_DATA_ROOTS=/home/zhaowl/Projects:/home/zhaowl/2026BatCN_NHZY_Yunnan
```

不要写 Windows 路径，也不要设置为 `/`。其余并发上限应按服务器资源调整，例如 96 线程服务器可保留默认值：

```bash
MAX_TOTAL_THREADS=96
MAX_QC_PARALLEL=8
MAX_ASSEMBLY_PARALLEL=4
MAX_THREADS_PER_MEGHIT=48
```

同时确认运行账号拥有输入目录的读取权限，以及 cleandata、assembly 的父目录写入权限：

```bash
ls -ld /home/zhaowl/2026BatCN_NHZY_Yunnan
```

## 4. 前台启动并测试

保持 `contig-ui` 已激活，在项目目录运行：

```bash
conda activate contig-ui
./.venv/bin/streamlit run app.py \
  --server.address 127.0.0.1 \
  --server.port 8501 \
  --server.headless true
```

另开一个 Linux 终端验证服务：

```bash
curl -I http://127.0.0.1:8501
```

正常时启动终端会显示本地访问地址，`curl` 会返回 HTTP 响应。前台测试通过后，按 `Ctrl+C` 停止，或继续保持运行。

## 5. Windows 通过 SSH 隧道访问

在 Windows PowerShell 执行（将 `fat1` 改为可解析的服务器地址或 IP）：

```powershell
ssh -N -L 8501:127.0.0.1:8501 zhaowl@fat1
```

不要关闭这个 PowerShell 窗口。随后在 Windows 浏览器打开：

```text
http://127.0.0.1:8501
```

界面中的路径始终是 Linux 服务器路径，例如 `/home/zhaowl/...`，不是 `C:\...`。

## 6. 常驻运行（可选，需 systemd 权限）

先确认 Conda 环境位置：

```bash
conda activate contig-ui
echo "$CONDA_PREFIX"
```

编辑 `deploy/contig-pipeline.service.example`，将 `User`、`WorkingDirectory`、`ExecStart` 改成你的实际路径；并添加 `PATH`，确保服务可找到 Conda 中的 fastp 和 MEGAHIT。例如 Conda 环境为 `/home/zhaowl/miniconda3/envs/contig-ui`：

```ini
User=zhaowl
WorkingDirectory=/home/zhaowl/Software/contig_pipeline
Environment=PATH=/home/zhaowl/miniconda3/envs/contig-ui/bin:/usr/local/bin:/usr/bin:/bin
Environment=LANG=C
Environment=LC_ALL=C
ExecStart=/home/zhaowl/Software/contig_pipeline/.venv/bin/streamlit run /home/zhaowl/Software/contig_pipeline/app.py --server.address 127.0.0.1 --server.port 8501 --server.headless true
```

由拥有 sudo 权限的管理员安装服务：

```bash
sudo cp deploy/contig-pipeline.service.example /etc/systemd/system/contig-pipeline.service
sudo systemctl daemon-reload
sudo systemctl enable --now contig-pipeline
sudo systemctl status contig-pipeline
```

查看运行日志：

```bash
sudo journalctl -u contig-pipeline -f
```

服务应始终以普通账号运行，不能以 root 运行。

## 7. 常见问题

### 仍显示 Python 3.6

说明没有激活新环境，或调用了错误解释器。依次执行：

```bash
conda activate contig-ui
which python
python --version
PYTHON_BIN=python bash install.sh
```

### 安装 pandas 时出现 `Meson`、`Cython` 或 GCC 4.8 编译错误

这表示 pip 没有使用二进制 wheel，而是在本地编译 pandas；不要在该服务器上为此升级 GCC。请先从本项目获取更新后的 `requirements.txt` 与 `install.sh`，然后在 `contig-ui` 环境内重新创建项目虚拟环境：

```bash
conda activate contig-ui
cd /home/zhaowl/Software/contig_pipeline
rm -rf .venv
PYTHON_BIN=python bash install.sh
```

若暂时不能更新项目文件，可在当前 `requirements.txt` 末尾增加一行 `pandas==2.2.3`，然后执行上述三条命令。pandas 2.2.3 为 CPython 3.10 提供 glibc 2.17+ 的 manylinux wheel。

### 安装 pyarrow 时出现 `libcst`、Rust compiler 或 `pyarrow-*.tar.gz`

这是同一类问题：pip 下载了 pyarrow 源码包，继而尝试构建其依赖。不要安装 Rust。更新项目的 `requirements.txt` 与 `install.sh` 后重建 `.venv`；新版会固定 `pyarrow==15.0.2` 并禁止所有源码构建：

```bash
conda activate contig-ui
cd /home/zhaowl/Software/contig_pipeline
rm -rf .venv
PYTHON_BIN=python bash install.sh
```

若暂时只能手工修改文件，在 `requirements.txt` 末尾追加：

```text
numpy==1.26.4
pyarrow==15.0.2
```

然后重新创建 `.venv`，并用 `./.venv/bin/python -m pip install --only-binary=:all: -r requirements.txt` 安装。pyarrow 15.0.2 提供 CPython 3.10、glibc 2.17+ 的 Linux wheel。

### `fastp` 或 `megahit` not in PATH

在启动网页服务之前激活 `contig-ui`。若 systemd 运行，检查 unit 中的 `Environment=PATH=.../envs/contig-ui/bin:...`。

### 网页打不开

先在 Linux 检查 `curl -I http://127.0.0.1:8501`，再确认 Windows 的 SSH 隧道仍在运行。第一版故意只监听 Linux 的 `127.0.0.1`，不要改为 `0.0.0.0`。

### 界面提示路径不在允许范围

修改 `config/pipeline.env` 的 `ALLOWED_DATA_ROOTS`，保存后重启 Streamlit。不要通过放宽到根目录 `/` 来绕过检查。
