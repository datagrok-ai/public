# Scripting changelog

## 1.1.4 (2026-11-04)

* Rebuilt the jkg_* kernel images with zero vulnerability findings: tornado 6.6 from PyPI in the gateway and Python environments (conda-forge stops at 6.5.10); the Octave kernel runs from the gateway's Python 3.11 instead of the Octave prefix's 3.10; the Node.js kernel moved from Node 18 (end of life) to Node 22, with a small zeromq 6 kernel replacing ijavascript's jmp/zeromq 5 transport (nel still runs the code, so scripts behave the same)

## 1.1.3 (2026-11-03)

* Rebuilt the jkg_* kernel images from master: jkg_octave clears its tornado HIGHs (tornado 6.5.10 from PyPI; conda-forge has no Python 3.10 build past 6.5.8), the other images pick up current base and dependency patches

## 1.1.2 (2026-10-14)

* [core#811](https://github.com/datagrok-ai/core/issues/811): Rebuilt the jkg_* kernel images from master: a script that a kernel does not start within 15 s is resent to a fresh kernel (Octave), instead of the worker reusing a wedged kernel until restarted

## 1.1.1 (2026-09-08)

* GROK-20743: Rebuilt jkg_nodejs — 1.1.0 shipped a js-api bundle without the Node runtime ports, so every Node.js script failed at bootstrap with "NoSuchMethodError: method not found: 'call'"; the other kernel images are the 1.1.0 builds under the new tag

## 1.1.0 (2026-08-12)

* GROK-18695: Rebuilt the jkg_* kernel images with remediated vulnerabilities — 357 findings (16 critical / 153 high) across the six images reduced to 5, none critical
* GROK-18695: Upgraded the Python science environment to Python 3.11 — mlflow 3.15, pyarrow 24, pytorch ≥2.6, transformers 5, pillow 12, pandas 2.3, numpy 1.26 (tensorflow 2.15 pip wheel keeps the Keras 2 API)
* GROK-18695: Replaced pycld3 (no builds past Python 3.8) with a cld3-compatible shim over conda-forge gcld3, keeping `import cld3` working
* GROK-18695: Removed unused packages from the Python environment (kyotocabinet, python-flair, jupyter/gunicorn/Flask leftovers, the django-pulling `image` package)
* GROK-18695: Modernized the Octave environment toolchain (python 3.10, current pip/setuptools) and dropped CPython's bundled ensurepip wheels from all conda environments
