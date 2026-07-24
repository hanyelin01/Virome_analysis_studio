#!/usr/bin/env bash
# Run on the Linux server from this directory: bash install.sh
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
command -v "$PYTHON_BIN" >/dev/null || { echo "Python executable not found: $PYTHON_BIN" >&2; exit 127; }
"$PYTHON_BIN" - <<'PY'
import sys
if sys.version_info < (3, 8):
    sys.stderr.write(
        "ERROR: This application requires Python 3.8 or newer; Python 3.10 is recommended.\n"
        f"Current interpreter: {sys.executable} ({sys.version.split()[0]})\n"
        "Create/activate a newer Conda environment, then run:\n"
        "  PYTHON_BIN=python bash install.sh\n"
    )
    raise SystemExit(3)
PY
[[ -f "$ROOT/config/pipeline.env" ]] || cp "$ROOT/config/pipeline.env.example" "$ROOT/config/pipeline.env"
if [[ -x "$ROOT/.venv/bin/python" ]]; then
  existing_version=$("$ROOT/.venv/bin/python" -c 'import sys; print("%s.%s" % sys.version_info[:2])')
  requested_version=$("$PYTHON_BIN" -c 'import sys; print("%s.%s" % sys.version_info[:2])')
  if [[ "$existing_version" != "$requested_version" ]]; then
    cat >&2 <<EOF
ERROR: Existing $ROOT/.venv uses Python $existing_version, but the selected interpreter uses Python $requested_version.
To avoid mixing packages from different Python versions, remove only this project's old virtual environment, then run this installer again:
  rm -rf "$ROOT/.venv"
EOF
    exit 3
  fi
fi
"$PYTHON_BIN" -m venv "$ROOT/.venv"
"$ROOT/.venv/bin/python" -m pip install --upgrade 'pip<25'
# Do not fall back to compiling any Python package from source on older compute
# servers. The requirements pin known compatible manylinux wheels.
"$ROOT/.venv/bin/python" -m pip install --only-binary=:all: -r "$ROOT/requirements.txt"
chmod +x "$ROOT/install.sh" "$ROOT/scripts/"*.sh
for tool in fastp megahit find realpath flock genomad checkv coverm vclust minimap2 samtools diamond taxonkit; do
  command -v "$tool" >/dev/null 2>&1 || echo "WARNING: $tool is not currently in PATH. The pipeline will not run until it is available."
done
echo "Installed. Edit $ROOT/config/pipeline.env, then start:"
echo "  $ROOT/.venv/bin/streamlit run $ROOT/app.py --server.address 127.0.0.1 --server.port 8501 --server.headless true"
