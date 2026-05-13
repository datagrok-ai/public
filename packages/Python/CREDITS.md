# Python — Third-Party Libraries

The `python` package is distributed under the MIT license that covers the rest
of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)). It
incorporates the open-source components listed below; this file reproduces the
attribution and notices required by their respective licenses.

The TypeScript bundle has no third-party runtime dependencies beyond the
platform-provided webpack externals. All third-party content is shipped via the
`dockerfiles/python` Docker image (Section 4) under permissive licenses.

---

## 1. Bundled in the published artifact (`dist/`)

None. All TypeScript runtime dependencies are either Datagrok-internal
(`@datagrok-libraries/*`, `datagrok-api`) and covered by the repo-wide MIT
`LICENSE.md`, or webpack externals provided by the platform host (Section 2).

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Python's `dist/` — they are provided once
by the platform host and shared across all packages. Listed here for full
transparency; their license texts ship with the platform itself.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

These images are built and distributed separately from the JavaScript bundle.
All listed components are under permissive licenses.

### `dockerfiles/python` — Python execution server

| Component    | Source                                     | License      |
|--------------|--------------------------------------------|--------------|
| Python 3.10  | https://www.python.org/                    | PSF License  |
| RDKit 2023.3.2 | https://www.rdkit.org/                   | BSD-3-Clause |
| pika         | https://github.com/pika/pika               | BSD-3-Clause |
| matplotlib   | https://matplotlib.org/                    | matplotlib (BSD-style) |
| NumPy        | https://numpy.org/                         | BSD-3-Clause |
| pandas       | https://pandas.pydata.org/                 | BSD-3-Clause |
| scikit-learn | https://scikit-learn.org/                  | BSD-3-Clause |
| websockets   | https://github.com/python-websockets/websockets | BSD-3-Clause |

The container also installs Debian system libraries (`libxrender1`,
`libxext6`, `libsm6`, `libglib2.0-0`, `libexpat1`) provided by the upstream
Debian distribution under their respective licenses.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack toolchain. These are
**not** redistributed as part of the published Python plugin and impose no
obligation on users of the plugin.
