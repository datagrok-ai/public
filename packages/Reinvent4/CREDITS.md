# Reinvent4 — Third-Party Libraries

The `@datagrok/reinvent4` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

The TypeScript bundle has only one third-party runtime JavaScript dependency
(`jszip`, MIT). The heavy lifting — REINVENT 4 molecular generation,
DockStream, AutoDock Vina docking — runs server-side inside the Docker image
(Section 4).

---

## 1. Bundled in the published artifact (`dist/`)

### JSZip (3.10.1)

ZIP file reader/writer used to package model / config folders that are
uploaded to the REINVENT 4 Docker container.

- Upstream: https://stuk.github.io/jszip/
- License: **MIT** (chosen from upstream's `MIT OR GPL-3.0-or-later` dual license)

```
Copyright (c) 2009-2016 Stuart Knightley, David Duponchel, Franz Buchinger, António Afonso

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Reinvent4's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

These images are built and distributed separately from the JavaScript bundle.
The container clones AstraZeneca's REINVENT 4 and DockStream repositories at
image-build time, installs them into a Python 3.10 virtualenv together with
their pinned dependencies (`requirements-macOS.lock`), downloads pretrained
model weights from a public AWS S3 bucket, and downloads the AutoDock Vina
binary distribution.

### `dockerfiles/reinvent` — REINVENT 4 / DockStream Flask service

| Component | Source | License |
|-----------|--------|---------|
| Python 3.10 | https://www.python.org/ | PSF License |
| **REINVENT 4** (commit `f814377`) | https://github.com/MolecularAI/REINVENT4 | **Apache-2.0** (AstraZeneca) |
| **DockStream** | https://github.com/MolecularAI/DockStream | **Apache-2.0** (AstraZeneca) |
| **AutoDock Vina 1.1.2** | https://vina.scripps.edu/ | **Apache-2.0** (The Scripps Research Institute) |
| Miniforge / conda-forge | https://github.com/conda-forge/miniforge | BSD-3-Clause |
| Flask | https://flask.palletsprojects.com/ | BSD-3-Clause |
| flask_cors | https://github.com/corydolphin/flask-cors | MIT |
| descriptastorus | https://github.com/bp-kelley/descriptastorus | BSD-3-Clause |
| toml | https://github.com/uiri/toml | MIT |
| pandas | https://pandas.pydata.org/ | BSD-3-Clause |
| pydantic 1.10.x | https://github.com/pydantic/pydantic | MIT |
| AWS CLI | https://github.com/aws/aws-cli | Apache-2.0 |

REINVENT 4 carries a NOTICE file
(<https://github.com/MolecularAI/REINVENT4/blob/main/LICENSE>) requiring
attribution to AstraZeneca; it is preserved within the cloned `REINVENT4`
checkout inside the container. Per Apache-2.0 §4, the small in-place
`sed`-based modifications applied to `Reinvent.py` and
`run_staged_learning.py` (replacing `model_dump()` with `model_config`) and
the DockStream `environment.yml` edits (Python version bump, pydantic pin)
are noted here.

The container also fetches Debian system packages (`libgomp1`, `swig`,
`cmake`, `build-essential`, `python3-dev`, `libeigen3-dev`, `zlib1g-dev`,
`libxrender1`, `libxext6`, `libsm6`, `libgl1`, etc.) from the upstream Debian
distribution under their respective licenses, and pretrained model
checkpoints from the public `s3://datagrok-data/models/reinvent` bucket.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack toolchain. These are
**not** redistributed as part of the published Reinvent4 plugin and impose no
obligation on users of the plugin.
