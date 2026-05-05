# TensorFlow.js — Third-Party Libraries

The `@datagrok/tensorflow-dev` package is distributed under the MIT license
that covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)). It
incorporates the open-source components listed below; this file reproduces the
attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (Apache-2.0, MIT). No copyleft (GPL/LGPL/MPL) component is
bundled into the published plugin (the GPL portion of jszip's dual MIT/GPLv3
license is **not** elected — see below).

---

## 1. Bundled in the published artifact (`dist/`)

### TensorFlow.js (`@tensorflow/tfjs` 4.x)

The full TensorFlow.js stack — `tfjs-core`, `tfjs-backend-cpu`,
`tfjs-backend-webgl`, `tfjs-converter`, `tfjs-data`, `tfjs-layers` — bundled
via the umbrella `@tensorflow/tfjs` npm package.

- Upstream: https://github.com/tensorflow/tfjs
- License: **Apache-2.0**

```
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

Copyright 2017–2024 Google LLC. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```

The full Apache-2.0 license text is available at the URL above. The
TensorFlow.js npm tarball does not ship a separate `NOTICE` file; the
copyright notice reproduced above appears in the per-file headers throughout
the upstream source. No modifications are made to TensorFlow.js by the
Datagrok plugin.

### JSZip (3.10.x)

ZIP file reader/writer used to package and unpackage TensorFlow.js model
artifacts (`model.json` + weight shards).

- Upstream: https://stuk.github.io/jszip/
- License: **MIT** (chosen from upstream's `MIT OR GPL-3.0-or-later` dual
  license)

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

### jszip-sync (3.2.x)

Fork of JSZip exposing a synchronous API; same upstream codebase as JSZip
above.

- Upstream: https://github.com/jhuckaby/jszip-sync
- License: **MIT** (chosen from upstream's `MIT OR GPL-3.0-or-later` dual
  license — same notice as JSZip; *Copyright (c) 2009-2016 Stuart Knightley,
  David Duponchel, Franz Buchinger, António Afonso*)

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into TensorFlow.js' `dist/` — they are
provided once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| OpenChemLib JS                   | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js        |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |

---

## 3. Fetched at runtime from third-party services

None.

---

## 4. Docker container images

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published TensorFlow.js plugin and impose no obligation on
users of the plugin.
