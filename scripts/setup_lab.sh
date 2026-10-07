#!/usr/bin/env bash
# Reproduce the dry lab exactly (Debian 12, Python 3.13, CPU only).
# Usage: bash scripts/setup_lab.sh
set -euo pipefail

cd "$(dirname "$0")/.."

# 1. Python environment
if ! command -v uv >/dev/null 2>&1; then
  python3 -m venv .venv
  ./.venv/bin/pip install -U pip
  ./.venv/bin/pip install -r scripts/requirements.txt
else
  uv venv --python 3.13 .venv
  uv pip install --python .venv/bin/python -r scripts/requirements.txt
fi

# 2. AutoDock Vina 1.2.5 static binary (no Boost build needed)
mkdir -p bin
if [ ! -x bin/vina ]; then
  curl -sL -o bin/vina \
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
  chmod +x bin/vina
fi

# 3. Report what we got
./.venv/bin/python - <<'PY'
import importlib
for m in ["rdkit", "numpy", "pandas", "scipy", "sklearn", "Bio", "openbabel", "meeko", "gemmi", "requests"]:
    mod = importlib.import_module(m)
    print(f"{m:10s} {getattr(mod, '__version__', '?')}")
PY
./bin/vina --version
echo "real CPU count: $(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)"
