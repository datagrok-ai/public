# CVMTests — Third-Party Libraries

The `@datagrok/cvm-tests` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

CVMTests is a service package containing automated CVM (Compute Virtual
Machine) tests. The bundled JavaScript dependencies are all under permissive
licenses (MIT, BSD-style). No copyleft component is bundled.

---

## 1. Bundled in the published artifact (`dist/`)

### ws (8.x)

WebSocket client/server used by the `grok_pipe` stress-test harness.

- Upstream: https://github.com/websockets/ws
- License: **MIT**

```
Copyright (c) 2011 Einar Otto Stangvik <einaros@gmail.com>
Copyright (c) 2013 Arnout Kazemier and contributors
Copyright (c) 2016 Luigi Pinca and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into CVMTests's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |

---

## 3. Docker containers (`dockerfiles/`)

Two service containers are shipped to exercise CVM container management. Both
images are built from `python:3.8-alpine3.18` and are published only as part
of test runs, not as a customer-facing product.

### `dockerfiles/cvmtests-docker-test1`

| Component        | License            |
|------------------|--------------------|
| Python 3.8       | PSF License        |
| Alpine Linux 3.18 base | MIT (busybox) + various permissive licenses for base utilities |
| Flask 2.3.3      | BSD-3-Clause       |

### `dockerfiles/cvmtests-docker-test2`

| Component        | License            |
|------------------|--------------------|
| Python 3.8       | PSF License        |
| Alpine Linux 3.18 base | MIT + various permissive licenses |
| Quart 0.18.4     | MIT                |
| hypercorn 0.14.4 | MIT                |
| websockets 11.0.3| BSD-3-Clause       |
| Werkzeug 2.3.7   | BSD-3-Clause       |

All listed components are under permissive licenses.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published CVMTests plugin and impose no obligation on users of
the plugin.
