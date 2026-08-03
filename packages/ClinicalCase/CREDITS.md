# ClinicalCase — Third-Party Libraries

The `@datagrok/clinical-case` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT, Apache-2.0). The Docker image referenced in
Section 3 fetches the CDISC Rules Engine `core` binary from upstream GitHub
releases — its license is reviewed below.

---

## 1. Bundled in the published artifact (`dist/`)

### jStat (1.9.x)

Statistical distributions library used in the SDTM/SEND analysis pipeline.

- Upstream: https://github.com/jstat/jstat
- License: **MIT**

```
Copyright (c) 2013 jStat

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

### Moment.js (2.29.x)

Date/time handling for clinical trial timelines.

- Upstream: https://momentjs.com/
- License: **MIT**

```
Copyright (c) JS Foundation and other contributors

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
```

### x2js (3.4.x)

XML ↔ JavaScript object conversion for SDTM XML payloads.

- Upstream: https://github.com/x2js/x2js
- Author: Axinom
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
> License for the specific language governing permissions and limitations
> under the License.

The full Apache-2.0 license text is available at the URL above and in
`node_modules/x2js/LICENSE`. The upstream tarball ships no `NOTICE` file.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into ClinicalCase's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |

---

## 3. Docker container (`dockerfiles/`) — CDISC Rules Engine service

The ClinicalCase package ships a single Docker image (`dockerfiles/Dockerfile`)
that runs the [CDISC Rules Engine](https://github.com/cdisc-org/cdisc-rules-engine)
`core` binary as a Celery worker. The image is built from `python:3.12-slim`
and pulls third-party software at build time:

| Component             | Source                                                                          | License                              |
|-----------------------|---------------------------------------------------------------------------------|--------------------------------------|
| CDISC Rules Engine    | https://github.com/cdisc-org/cdisc-rules-engine (latest GitHub release `core`)  | MIT (CDISC, Inc.)                    |
| Python 3.12 base      | Debian-slim base image, https://www.python.org/                                 | PSF License (Python) + Debian (base) |
| Celery                | https://github.com/celery/celery                                                | BSD-3-Clause                         |
| `datagrok-celery-task`| Internal Datagrok package (PyPI)                                                | MIT (covered by repo umbrella)       |
| `datagrok-api` (py)   | Internal Datagrok package (PyPI)                                                | MIT (covered by repo umbrella)       |
| filelock              | https://github.com/tox-dev/filelock                                             | Unlicense                            |
| jq, curl, unzip       | Debian apt packages                                                             | MIT / MIT / GPLv3 respectively       |

`unzip` is licensed under the Info-ZIP license (BSD-style) and `jq` under MIT.
The container's installed `unzip` is taken from Debian and is only used at
build time to extract the CDISC Rules Engine release archive — it is included
in the running image but is not invoked at runtime.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published ClinicalCase plugin and impose no obligation on users
of the plugin.
