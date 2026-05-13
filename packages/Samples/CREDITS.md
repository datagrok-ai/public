# Samples — Third-Party Libraries

The `@datagrok/samples` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This package contains sample data connectors, scripts, and notebooks. It does
not bundle any third-party JavaScript runtime libraries beyond what the
Datagrok platform already provides as webpack externals. The package does
however ship one Docker container that packages a small Python ML demo (see
Section 4).

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either listed in Section 2 (provided by the
platform) or are Datagrok-internal libraries covered by the repo-wide
`LICENSE.md`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Samples' `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License    | Upstream                          |
|-----------|---------|------------|-----------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs |

---

## 3. Fetched at runtime from third-party services

None.

---

## 4. Docker container (`dockerfiles/`) — Python KNN demo worker

The Samples package ships a single Docker image (`dockerfiles/Dockerfile`)
that runs a Celery worker exposing a small Python K-nearest-neighbours model
as a Datagrok function. The image is built from `mambaorg/micromamba:1.5.3`
and pulls third-party software at build time:

| Component                    | Source                                             | License                |
|------------------------------|----------------------------------------------------|------------------------|
| micromamba (mamba)           | https://github.com/mamba-org/mamba                | BSD-3-Clause           |
| Python 3.8                   | https://www.python.org/                            | PSF License            |
| NumPy                        | https://numpy.org/                                  | BSD-3-Clause           |
| pandas                       | https://pandas.pydata.org/                          | BSD-3-Clause           |
| scikit-learn                 | https://scikit-learn.org/                           | BSD-3-Clause           |
| Celery                       | https://docs.celeryq.dev/                           | BSD-3-Clause           |
| datagrok-celery-task         | https://pypi.org/project/datagrok-celery-task/      | MIT (Datagrok)         |
| Debian base + build-essential| https://www.debian.org/                             | various permissive (GPL/LGPL/BSD/MIT) |

Standard Debian build tooling (`gcc`, `g++`, `gfortran`, `build-essential`) is
installed at image-build time; per the GNU GCC Runtime Library Exception,
binaries produced by GCC are not subject to the GPL.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Samples plugin and impose no obligation on users of
the plugin.

The dev dependencies on other Datagrok plugins (`@datagrok/api-tests`,
`@datagrok/dendrogram`) are MIT — covered by the repo-wide `LICENSE.md`.
