# Flow — Third-Party Libraries

The `@datagrok/flow` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is bundled
into the published plugin artifact.

> **Note on this file.** Flow is built without `node_modules` checked in
> locally; the per-package LICENSE files were not available at the time this
> credits file was generated. The license texts below are reproduced from the
> canonical upstream sources for the version pinned in `package.json`. The
> SPDX identifiers were taken from the same `package.json` declarations.

---

## 1. Bundled in the published artifact (`dist/`)

### litegraph.js (0.7.x)

Node-based graph editor used as the canvas back-end for the visual function-flow
designer.

- Upstream: https://github.com/jagenjo/litegraph.js
- License: **MIT** (per upstream `LICENSE`)

```
The MIT License (MIT)

Copyright (c) 2013 Javi Agenjo (javi.agenjo@gmail.com)

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

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Flow's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Flow plugin and impose no obligation on users of the
plugin.
