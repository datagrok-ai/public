# Dendrogram — Third-Party Libraries

The `@datagrok/dendrogram` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is
bundled into the published Dendrogram plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### fastest-levenshtein (1.0.x)

Levenshtein-distance implementation used for sequence-distance computation
inside the dendrogram clustering routines.

- Upstream: https://github.com/ka-weihe/fastest-levenshtein
- License: **MIT**

```
MIT License

Copyright (c) 2020 Kasper Unn Weihe

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

(`@webgpu/types`, `@types/wu`, and `@types/node` are TypeScript-only types
stubs — no runtime code is shipped from them. `file-loader`, `worker-loader`,
and `source-map-loader` are webpack build plugins and produce no bundled
runtime code.)

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Dendrogram's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| OpenChemLib JS                   | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js        |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |
| ExcelJS                          | 4.x     | MIT          | https://github.com/exceljs/exceljs                |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Dendrogram plugin and impose no obligation on users
of the plugin.
