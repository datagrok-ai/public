# GenomeBrowser — Third-Party Libraries

The `@datagrok/genome-browser` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, Apache-2.0). No copyleft (GPL/LGPL/MPL) component
is bundled into the published GenomeBrowser plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### JBrowse 2 — @jbrowse/react-linear-genome-view (3.x)

Embeddable linear-genome-view React component from the JBrowse 2 project.
This is the core viewer rendered inside the GenomeBrowser plugin and pulls
in many `@jbrowse/*` subpackages (core, plugin-* modules, etc.) — each of
which ships its own MIT or Apache-2.0 LICENSE under
`node_modules/@jbrowse/*/`.

- Upstream: https://jbrowse.org/ — https://github.com/GMOD/jbrowse-components
- License: **Apache-2.0** (and **MIT** for some subpackages — see individual
  `node_modules/@jbrowse/*/LICENSE` files; `react-linear-genome-view`'s own
  `package.json` declares MIT)

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

JBrowse 2 transitively bundles a large number of bioinformatics file-format
parsers (BAM, CRAM, VCF, BigWig, GFF3, etc.), Material UI (MIT), MobX (MIT),
mobx-state-tree (MIT) and many others — each carrying its own permissive
license under `node_modules/`.

### React (19.x) and React DOM (19.x)

Required peer dependencies of the JBrowse component.

- Upstream: https://react.dev/ — https://github.com/facebook/react
- License: **MIT**

```
MIT License

Copyright (c) Meta Platforms, Inc. and affiliates.

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

### buffer (6.0.x)

Browser-side `Buffer` polyfill needed by JBrowse's binary file-format parsers.

- Upstream: https://github.com/feross/buffer
- License: **MIT**

```
The MIT License (MIT)

Copyright (c) Feross Aboukhadijeh, and other contributors.

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

These libraries are not bundled into GenomeBrowser's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| RxJS      | 7.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

GenomeBrowser overrides `rxjs` to 7.x (vs the platform default 6.x) via the
`overrides` field in `package.json`.

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

JBrowse 2 viewers can be configured to fetch genome assemblies, tracks (BAM,
CRAM, VCF, BigWig, GFF3), reference sequences, and ontology terms from
third-party servers (UCSC, Ensembl, NCBI, etc.) at runtime, depending on the
user's session configuration. These resources are **not** bundled into
GenomeBrowser's `dist/`. Each provider's data-use terms apply.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published GenomeBrowser plugin and impose no obligation on
users of the plugin.

Babel, Vite, ESLint, fork-ts-checker-webpack-plugin and similar tooling are
build-time only (all MIT or BSD) and not included in the runtime artifact.
