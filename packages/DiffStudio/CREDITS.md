# Diff Studio — Third-Party Libraries

The `@datagrok/diff-studio` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is
bundled into the published Diff Studio plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### CodeMirror 6 (`codemirror`, `@codemirror/autocomplete`, `@codemirror/lang-markdown`, `@codemirror/lang-python`, `@codemirror/language-data`, `@codemirror/state`)

Modular code editor framework used as the IVP-formula editor inside the Diff
Studio app. Each `@codemirror/*` package and the umbrella `codemirror` 6.x
package ships with its own copy of the same MIT notice, reproduced below.

- Upstream: https://codemirror.net/
- License: **MIT**

```
MIT License

Copyright (C) 2018-2021 by Marijn Haverbeke <marijn@haverbeke.berlin> and others

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

CodeMirror 6 is composed of many sibling packages; the plugin pulls in the
six packages listed above directly and additional `@codemirror/*` and
`@lezer/*` packages transitively. All of them ship under the same MIT notice;
see `node_modules/@codemirror/*/LICENSE` and `node_modules/@lezer/*/LICENSE`
for per-package copies.

### diff-grok (1.2.0)

Zero-dependency TypeScript ODE solver library implementing CVODE, LSODA,
Rosenbrock-Wanner and Runge-Kutta methods. This is the numerical engine
underlying Diff Studio.

- Upstream: https://www.npmjs.com/package/diff-grok
- License: **MIT**

```
MIT License

Copyright (c) 2026 Datagrok

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

The `diff-grok` package itself includes a TypeScript port of SUNDIALS CVODE
v7.5.0; the upstream SUNDIALS suite (Hindmarsh et al., LLNL) is distributed
under BSD-3-Clause. See `node_modules/diff-grok/THIRD_PARTY_LICENSES` for the
full upstream attribution preserved by the diff-grok authors.

### @babel/runtime (7.27.x)

Babel runtime helpers (regenerator runtime, helper functions) injected by the
build to support ES2018+ down-leveled code paths.

- Upstream: https://babeljs.io/
- License: **MIT**

```
MIT License

Copyright (c) 2014-present Sebastian McKenzie and other contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

(`css-loader` and `style-loader` are listed in `dependencies` but are
build-time-only webpack plugins and produce no bundled runtime code.)

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Diff Studio's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| datagrok-api                     | 1.x     | MIT          | https://github.com/datagrok-ai/public             |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Diff Studio plugin and impose no obligation on users
of the plugin.
