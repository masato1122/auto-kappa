
python setup.py sdist
# Install auto-kappa and its dependencies in a virtual environment
# Install matbench-discovery
git clone https://github.com/janosh/matbench-discovery --depth 1
uv pip install -e ./matbench-discovery
# Install MLIPs
# fairchem-core-1.10.0
uv pip install -e packages/fairchem-core
# Install auto-kappa
uv pip install dist/auto_kappa-1.1.0.tar.gz

rm -r dist

