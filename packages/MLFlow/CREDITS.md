# MLflow — Third-Party Libraries

The `@datagrok/mlflow` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

The TypeScript plugin itself bundles no third-party JavaScript runtime
dependencies — it only uses the platform's webpack externals. The package's
`dockerfiles/` directory ships a Docker image that hosts an MLflow tracking
server; third-party content baked into that image is itemized below.

---

## 1. Bundled in the published artifact (`dist/`)

None. MLflow's TypeScript plugin code is the only application-level code in
the bundle (apart from the webpack runtime). All external libraries it uses
are platform externals listed in Section 2.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into MLflow's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Docker container images (`dockerfiles/`)

### `dockerfiles/mlflow` — MLflow tracking server

The image is built `FROM ghcr.io/mlflow/mlflow:latest` and adds standard
build tooling, `pyenv`, `boto3`, and additional Python interpreters. All
listed components are under permissive licenses.

| Component                            | Source                                     | License      |
|--------------------------------------|--------------------------------------------|--------------|
| MLflow                               | https://github.com/mlflow/mlflow           | Apache-2.0   |
| Python (3.10 / 3.11 / 3.12 via pyenv)| https://www.python.org/                    | PSF          |
| pyenv                                | https://github.com/pyenv/pyenv             | MIT          |
| boto3 / botocore                     | https://github.com/boto/boto3              | Apache-2.0   |
| virtualenv                           | https://virtualenv.pypa.io/                | MIT          |
| Debian base packages (curl, git,     | https://www.debian.org/                    | various      |
|   gcc, make, libssl-dev, …)          |                                            | (mostly      |
|                                      |                                            | GPL/LGPL)    |
| MeCab IPA dictionary (UTF-8)         | https://taku910.github.io/mecab/           | BSD-3-Clause |

The Debian base image and the MLflow upstream image both include GNU/Linux
system libraries under copyleft licenses (glibc — LGPL; coreutils, bash,
apt — GPL). These are standard OS components used as-is and are subject to
the obligations of their upstream licenses if the image itself is
redistributed.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published MLflow plugin and impose no obligation on users of
the plugin.
