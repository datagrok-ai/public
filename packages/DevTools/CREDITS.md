# DevTools — Third-Party Libraries

The `@datagrok/dev-tools` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

DevTools provides Test Manager and Dev Panel for plugin developers. No
copyleft component is bundled.

---

## 1. Bundled in the published artifact (`dist/` and `common/`)

### CodeMirror (5.65.x)

Source code editor component used in the Dev Panel for inspecting / editing
function bodies. CodeMirror is **not** bundled into Webpack output (it is
declared as a webpack external) but the package's `sources` field ships
`common/codemirror/codemirror.js` and `common/codemirror/codemirror.css`
alongside `dist/`, so the platform serves them at runtime.

- Upstream: https://codemirror.net/5/
- License: **MIT**

```
MIT License

Copyright (C) 2017 by Marijn Haverbeke <marijn@haverbeke.berlin> and others

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

These libraries are not bundled into DevTools's webpack output — they are
provided once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| CodeMirror 5                     | 5.65.x  | MIT          | https://codemirror.net/5/                         |

CodeMirror appears both here (as a webpack external resolved against the global
`CodeMirror` object) and in Section 1 (the actual JS/CSS files shipped via the
`sources` field). Both rows describe the same upstream package; only the
distribution mechanism differs.

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published DevTools plugin and impose no obligation on users of
the plugin.
