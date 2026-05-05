# SQLite — Third-Party Libraries

The `@datagrok/sqlite` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)). It
incorporates the open-source components listed below; this file reproduces the
attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, public-domain). No copyleft (GPL/LGPL/MPL) component
is bundled.

---

## 1. Bundled in the published artifact (`dist/`)

### sql.js (SQLite compiled to WebAssembly)

Used to read and preview SQLite database files entirely in the browser. The
two artifacts checked into `src/` (`sql-wasm.js`, `sql-wasm.wasm`) are copied
into `dist/` at build time.

- Upstream: https://github.com/sql-js/sql.js
- License: **MIT** (sql.js wrapper) **AND public-domain** (the SQLite C
  source compiled into the bundled WASM module).

The sql.js npm distribution does not include a per-file copyright header in
the JavaScript glue (the file `src/sql-wasm.js` begins with implementation
comments only). The upstream `sql.js` LICENSE is reproduced verbatim from the
project repository:

```
The MIT License (MIT)

Copyright (c) 2017 Ophir LOJKINE

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

The SQLite engine itself, compiled into `sql-wasm.wasm`, is in the **public
domain** (https://www.sqlite.org/copyright.html). The Emscripten runtime
embedded in `sql-wasm.js` is **MIT** (Copyright (c) Mozilla Foundation,
Emscripten authors).

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into SQLite's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| (none)    | —       | —       | (the package only externalizes `datagrok-api/{dg,grok,ui}`) |

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
as part of the published SQLite plugin and impose no obligation on users of
the plugin.

The `@types/sql.js` dev dependency contains TypeScript type declarations only
(MIT, DefinitelyTyped contributors); it ships no executable code into `dist/`.
