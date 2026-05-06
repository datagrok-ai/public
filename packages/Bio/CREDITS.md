# Bio — Third-Party Libraries

The `@datagrok/bio` package is distributed under the MIT license that covers the
rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)). It
incorporates the open-source components listed below; this file reproduces the
attribution and notices required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT, Apache-2.0, BSD-3-Clause). No copyleft
(GPL/LGPL/MPL) component is bundled into the published Bio plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### immunum (1.1.0)

ENPICOM's antibody/TCR numbering library, compiled to WebAssembly. The
WASM binary is bundled into `dist/`; the wasm-bindgen JS glue is reused from
[`src/utils/antibody-numbering/immunum-glue.js`](src/utils/antibody-numbering/immunum-glue.js),
which is a port of the upstream `node_modules/immunum/immunum.js` modified
only to replace the Node-only `require('fs').readFileSync` top-level loader
with an explicit `initImmunum(bytes)` entry point so it works in browsers and
web workers. The modification is noted in the file header per MIT terms.

- Upstream: https://github.com/ENPICOM/immunum — https://immunum.enpicom.com
- License: **MIT**

```
MIT License

Copyright (c) 2026 ENPICOM

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

### @biowasm/aioli (3.x)

WebAssembly runtime that drives multiple sequence alignment (kalign). The
aioli library itself is bundled; the kalign tool it loads is fetched at
runtime — see Section 3 below.

- Upstream: https://github.com/biowasm/aioli
- License: **MIT**

```
MIT License

Copyright (c) 2018 Robert Aboukhalil

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

### Ajv (8.x) and ajv-errors (3.x)

JSON-schema validator used for monomer-library JSON validation.

- Upstream: https://ajv.js.org/ — https://github.com/ajv-validator/ajv-errors
- License: **MIT**

```
The MIT License (MIT)

Copyright (c) 2015-2021 Evgeny Poberezkin

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

(`ajv-errors` carries an analogous MIT notice — *Copyright (c) 2017 Evgeny
Poberezkin*.)

### fastest-levenshtein (1.0.x)

Levenshtein-distance implementation used in sequence comparison code paths.

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

### umap-js (1.4.x)

UMAP dimensionality reduction used by sequence-space and similarity viewers.
The npm package metadata declares MIT but the LICENSE file shipped in the
tarball is **Apache-2.0**, so we comply with the more conservative Apache-2.0
terms (preserve copyright/attribution; mark modifications, of which there are
none).

- Upstream: https://github.com/PAIR-code/umap-js
- Author: Andy Coenen (Google PAIR)
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
> License for the specific language governing permissions and limitations
> under the License.

The full Apache-2.0 license text is available at the URL above and in
`node_modules/umap-js/LICENSE`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Bio's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| OpenChemLib JS                   | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js        |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |

OpenChemLib JS is also referenced by Bio's `package.json` `sources`. The
upstream `openchemlib-js` npm package is distributed under BSD-3-Clause:
*Copyright (c) 2015-2017, cheminfo.*

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

### kalign (via @biowasm/aioli)

`@biowasm/aioli` loads the **kalign** multiple-sequence-alignment tool from
the public biowasm CDN (`https://biowasm.com/cdn/v3/...`) on demand. kalign
itself is **not** bundled into Bio's `dist/` and is **not** redistributed by
the Datagrok platform. Users that load it transitively are subject to
kalign's upstream license (GPLv3).

- Upstream: https://github.com/TimoLassmann/kalign
- License: GPLv3 (used as a runtime-fetched, non-bundled dependency)

---

## 4. Docker container (`dockerfiles/`) — PepSeA service

The Bio package ships a single Docker image (`dockerfiles/Dockerfile`) that
provides a PepSeA-based MSA service for HELM peptide sequences. The image is
built from `datagrok/python` and pulls third-party software at build time:

| Component     | Source                                          | License                |
|---------------|-------------------------------------------------|------------------------|
| PepSeA        | https://github.com/Merck/PepSeA (main branch)   | MIT (Merck & Co.)      |
| MAFFT         | https://mafft.cbrc.jp/ (`mafft_7.520-1` .deb)   | MAFFT license (BSD-style, free for any use with attribution) |
| FastAPI       | https://github.com/tiangolo/fastapi             | MIT                    |
| uvicorn       | https://www.uvicorn.org/                        | BSD-3-Clause           |
| ujson         | https://github.com/ultrajson/ultrajson          | BSD-3-Clause           |
| Python stdlib | https://www.python.org/                         | PSF License            |

The `Dockerfile` patches a few lines of PepSeA's `api.py` at build time to
add a `/distout` endpoint and JSON error-handling middleware; per MIT terms
those modifications are noted here.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published Bio plugin and impose no obligation on
users of the plugin.

The peer/devDependencies on other Datagrok plugins (`@datagrok/chem`,
`@datagrok/dendrogram`, `@datagrok/eda`, `@datagrok/helm`,
`@datagrok/peptides`) are MIT — covered by the repo-wide `LICENSE.md`. Each
of those plugins maintains its own `CREDITS.md` (or will) for its own
third-party content.
