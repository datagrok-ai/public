# PowerPack — Third-Party Libraries

The `@datagrok/power-pack` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, Apache-2.0). No copyleft (GPL/LGPL/MPL) component
is bundled into the published PowerPack plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### CodeMirror 6 (`codemirror`, `@codemirror/autocomplete`, `@codemirror/lang-javascript`, `@codemirror/language`, `@codemirror/search`, `@codemirror/state`, `@codemirror/view`)

Code editor framework used by PowerPack's expression-editor / formula-lines
input.

- Upstream: https://codemirror.net/ — https://github.com/codemirror
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

### core-js (3.x)

ECMAScript polyfill suite (Promise, Symbol, iterators, etc.).

- Upstream: https://github.com/zloirock/core-js
- License: **MIT**

```
Copyright (c) 2013–2025 Denis Pushkarev (zloirock.ru)
Copyright (c) 2025–2026 CoreJS Company (core-js.io)

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

### ExcelJS (4.3.x)

Excel (`.xlsx`) workbook reader/writer used for the Excel-import path. Run
inside a web worker (`src/workers/exceljs-worker.ts`).

- Upstream: https://github.com/exceljs/exceljs
- License: **MIT**

```
The MIT License (MIT)

Copyright (c) 2014-2019 Guyon Roche

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

### Fuse.js (6.x)

Lightweight fuzzy-search library used by the spotlight / search-template
matcher.

- Upstream: https://github.com/krisk/Fuse
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy
> of the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

> Copyright 2017 Kirollos Risk

The full Apache-2.0 license text is available at the URL above and in
`node_modules/fuse.js/LICENSE`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into PowerPack's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

None.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI, the TypeScript / webpack / Babel toolchain, and the
`yaml` parser used at build time. These are **not** redistributed as part of
the published PowerPack plugin and impose no obligation on users of the
plugin.
