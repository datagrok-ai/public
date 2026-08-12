# Scripting changelog

## 1.1.0 (2026-08-12)

* GROK-18695: Rebuilt the jkg_* kernel images with remediated vulnerabilities — 357 findings (16 critical / 153 high) across the six images reduced to 5, none critical
* GROK-18695: Upgraded the Python science environment to Python 3.11 — mlflow 3.15, pyarrow 24, pytorch ≥2.6, transformers 5, pillow 12, pandas 2.3, numpy 1.26 (tensorflow 2.15 pip wheel keeps the Keras 2 API)
* GROK-18695: Replaced pycld3 (no builds past Python 3.8) with a cld3-compatible shim over conda-forge gcld3, keeping `import cld3` working
* GROK-18695: Removed unused packages from the Python environment (kyotocabinet, python-flair, jupyter/gunicorn/Flask leftovers, the django-pulling `image` package)
* GROK-18695: Modernized the Octave environment toolchain (python 3.10, current pip/setuptools) and dropped CPython's bundled ensurepip wheels from all conda environments
