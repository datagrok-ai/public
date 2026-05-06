# Docking — Third-Party Libraries

The `@datagrok/docking` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

The JavaScript portion of the plugin bundles **no third-party runtime
JavaScript** into the published artifact: every runtime dependency is either
internal (Datagrok platform / `@datagrok-libraries/*`) or provided by the
platform host as a webpack external.

The plugin's heavy lifting is performed inside a Docker container that ships
AutoDock Suite, AutoDock-GPU, and AutoDockTools. Those components carry their
own licenses, reviewed in Section 3 below. **The Docker image fetches GPL
software (AutoDockTools, AutoDock 4 / Vina, AutoDock-GPU) from upstream
sources at build time** — operators publishing the image must comply with
those upstream terms when redistributing.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Docking's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| OpenChemLib JS                   | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js        |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |
| NGL Viewer                       | 2.x     | MIT          | https://github.com/arose/ngl                      |

NGL Viewer is exposed as the `NGL` global by the platform host (`@datagrok/biostructure-viewer`) and is referenced as a webpack external from this plugin.

---

## 2. Docker container (`dockerfiles/autodock`) — AutoDock service

The Docking package ships a single Docker image (`dockerfiles/autodock/Dockerfile`)
that runs an AutoDock 4 / AutoDock-GPU docking service behind a Flask HTTP
endpoint. The image is built from `ubuntu:20.04` and pulls third-party
software at build time:

| Component                  | Source                                                                                | License                                          |
|----------------------------|---------------------------------------------------------------------------------------|--------------------------------------------------|
| Ubuntu 20.04 base          | https://hub.docker.com/_/ubuntu                                                       | Various permissive (apt packages)                |
| Miniforge / conda          | https://github.com/conda-forge/miniforge                                              | BSD-3-Clause                                     |
| mamba                      | https://github.com/mamba-org/mamba                                                    | BSD-3-Clause                                     |
| Python 2.7                 | https://www.python.org/                                                               | PSF License                                      |
| autodocktools-prepare      | https://github.com/insilichem/autodocktools_prep                                      | LGPL-2.1 (insilichem fork) wrapping AutoDockTools (Mozilla MPL-1.1 / LGPL by CCSB Scripps) |
| AutoDock Suite 4.2.6       | https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar | **GPLv2** (CCSB, Scripps Research Institute)     |
| AutoDock-GPU 1.5.3         | https://github.com/ccsb-scripps/AutoDock-GPU/releases/tag/v1.5.3                      | **GPLv2** (CCSB, Scripps Research Institute)     |
| `prepare_dpf4.py`          | https://github.com/MolecularFlipbook/FlipbookApp (mirrored AutoDockTools utility)     | MPL-1.1 / LGPL (AutoDockTools)                   |
| Flask                      | https://github.com/pallets/flask                                                      | BSD-3-Clause                                     |
| futures (Python)           | https://pypi.org/project/futures/                                                     | PSF License                                      |
| awscli (build-only)        | https://github.com/aws/aws-cli                                                        | Apache-2.0                                       |
| ocl-icd, opencl-headers, libgomp1, clinfo | Debian apt packages                                                  | LGPL / Apache-2.0 / MIT (mixed)                  |
| Datagrok-curated docking inputs | `s3://datagrok-data/autodock` (copied into `/app` at build time, Datagrok-owned) | MIT (covered by repo umbrella)                   |

### GPL note

AutoDock 4, AutoDock Vina (if used) and AutoDock-GPU are distributed by the
Center for Computational Structural Biology (CCSB) at The Scripps Research
Institute under the **GNU General Public License v2**. The Docker image
fetches their published binaries; redistributing this image (e.g. publishing
it to a public registry) requires complying with the GPLv2 — including
making the corresponding source available on request.

The autodocktools-prepare conda package and the mirrored
`prepare_dpf4.py` utility carry the upstream AutoDockTools license terms
(Mozilla Public License 1.1 / LGPL); these are **not** restated here, but
are preserved inside the conda environment and in the script's own header.

`autodock.py` (in this repository) is the Datagrok-authored wrapper that
exposes a Flask HTTP API in front of the GPL binaries; it is MIT-licensed
under the repo umbrella and links to AutoDock as a separate process, not as
linked code.

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Docking plugin and impose no obligation on users of
the plugin.

The peer/devDependency on `@datagrok/biostructure-viewer` is MIT — covered by
the repo-wide `LICENSE.md`. That plugin maintains its own `CREDITS.md` for its
own third-party content (NGL Viewer and friends).
