# Pyodide — Third-Party Libraries

The `@datagrok/pyodide` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT). Pyodide itself (Mozilla Public License 2.0) is
**not** bundled — it is loaded at runtime from a public CDN by the package's
web worker; see Section 3.

---

## 1. Bundled in the published artifact (`dist/`)

### uuid (10.x)

RFC-compliant UUID generation, used to correlate worker requests/responses
between the main thread and the Pyodide web worker.

- Upstream: https://github.com/uuidjs/uuid
- License: **MIT**

```
The MIT License (MIT)

Copyright (c) 2010-2020 Robert Kieffer and other contributors

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Pyodide's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

### Pyodide (0.27.7)

The package's web worker (`src/worker.js`) loads Pyodide and its CPython
WebAssembly runtime from the public jsDelivr CDN at runtime via
`importScripts("https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.js")`.
Pyodide is **not** bundled into Pyodide's `dist/` and is **not**
redistributed by the Datagrok platform; users that load it transitively are
subject to Pyodide's upstream license.

- Upstream: https://pyodide.org/ — https://github.com/pyodide/pyodide
- License: **Mozilla Public License 2.0** (file-level weak copyleft)

The MPL-2.0 license text is available at
https://www.mozilla.org/en-US/MPL/2.0/ and in the upstream Pyodide
distribution. CPython itself is bundled within Pyodide under the **PSF
License**, and the various scientific Python packages loaded on demand
(NumPy, pandas, etc.) ship under their own permissive licenses
(BSD-3-Clause, MIT) inside the Pyodide distribution.

### uuidv4 helper

The worker also loads `https://cdnjs.cloudflare.com/ajax/libs/uuid/8.2.0/uuidv4.min.js`
from cdnjs at runtime — the same `uuid` library credited in Section 1
(MIT, *Copyright (c) 2010-2020 Robert Kieffer and other contributors*).

---

## 4. Docker container images (`dockerfiles/`)

None. The Pyodide package runs entirely in the browser.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack toolchain. These are
**not** redistributed as part of the published Pyodide plugin and impose no
obligation on users of the plugin.
