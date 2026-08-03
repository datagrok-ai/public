# Mol Track — Third-Party Libraries

The `@datagrok/moltrack` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) JavaScript
component is bundled into the published Mol Track plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### randexp (0.5.x)

Generates a random string that matches a given regular expression — used
to materialize sample identifiers from regex templates.

- Upstream: http://fent.github.io/randexp.js/
- License: **MIT**

```
MIT License

Copyright (C) 2011 by fent

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

These libraries are not bundled into Mol Track's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

---

## 3. Docker container images (`dockerfiles/`)

### `dockerfiles/` — Mol Track registration service

The image is built `FROM python:3.11-slim` and installs the
[`dg-mol-track`](https://pypi.org/project/dg-mol-track/) Python package, which
provides the chemistry-registration backend.

| Component                       | Source                                         | License                  |
|---------------------------------|------------------------------------------------|--------------------------|
| dg-mol-track                    | https://pypi.org/project/dg-mol-track/         | MIT (Datagrok-published) |
| RDKit (transitive via dg-mol-track) | https://www.rdkit.org/                     | BSD-3-Clause             |
| FastAPI / Starlette             | https://fastapi.tiangolo.com/                  | MIT                      |
| SQLAlchemy                      | https://www.sqlalchemy.org/                    | MIT                      |
| psycopg / asyncpg               | https://www.psycopg.org/ / https://magicstack.github.io/asyncpg/ | LGPL-3.0+ / Apache-2.0 |
| pydantic                        | https://github.com/pydantic/pydantic           | MIT                      |
| uvicorn                         | https://www.uvicorn.org/                       | BSD-3-Clause             |
| jq (CLI)                        | https://stedolan.github.io/jq/                 | MIT                      |
| Python stdlib                   | https://www.python.org/                        | PSF                      |
| Debian slim base packages       | https://www.debian.org/                        | various (mostly GPL/LGPL)|

The Debian base image includes GNU/Linux system libraries under copyleft
licenses (glibc — LGPL; coreutils, bash, apt — GPL). These are standard OS
components used as-is and are subject to the obligations of their upstream
licenses if the image itself is redistributed. `psycopg2` (LGPL-3.0+) is
linked dynamically and used as a library — LGPL obligations apply if the
image is redistributed.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Mol Track plugin and impose no obligation on users of
the plugin.
