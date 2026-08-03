# EDA — Third-Party Libraries

The `@datagrok/eda` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, Apache-2.0, BSD-3-Clause). No copyleft (GPL/LGPL/MPL)
component is bundled into the published EDA plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### XGBoost (linked into `wasm/XGBoostAPI.wasm`)

Gradient-boosting library compiled to WebAssembly. The C++ XGBoost source
(upstream tag **v3.3.0**, with a minimal-build patch from
`wasm/xgboost/patches/` that trims unused components) is linked into
`wasm/XGBoostAPI.wasm` at build time (Emscripten, see `wasm/xgboost/build.ps1`;
XGBoost sources come from an external clone); the resulting
`wasm/XGBoostAPI.js` loader + `.wasm` binary are bundled by webpack into the
published artifact (the binary is emitted to `dist/`). XGBoost bundles
dmlc-core (also Apache-2.0, same copyright holders), which is linked in as
part of the build.

- Upstream: https://xgboost.ai/ — https://github.com/dmlc/xgboost
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

The full Apache-2.0 license text is reproduced in `THIRD_PARTY_LICENSES.txt`
and available at the URL above. Copyright is held by *XGBoost contributors*
(`Copyright (c) 2019 by Contributors`); upstream ships no NOTICE file.

### sci-comp-ml (`wasm/sci_comp_ml_bg.wasm`)

In-house Rust implementations of PCA (NIPALS), PLS1, softmax classification
and linear regression / Elastic Net, compiled to WebAssembly with
`wasm-pack`. The source lives in the separate `sci-comp-rust` repository and
is covered by the repo-wide MIT umbrella; it is listed here because the
compiled artifact (`wasm/sci_comp_ml.js` glue + `wasm/sci_comp_ml_bg.wasm`)
is bundled. The numerical kernels use no external crates; the WASM-boundary
dependencies (`wasm-bindgen`, `serde`, `serde-wasm-bindgen`, `js-sys`) are
**MIT** or **MIT OR Apache-2.0** (MIT elected).

### @keckelt/tsne (1.0.x)

t-SNE implementation in TypeScript, used by the dimensionality-reduction tools.

- Upstream: https://github.com/keckelt/tsne
- License: **MIT**

```
MIT License

Copyright (c) 2021 Klaus Eckelt

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

UMAP dimensionality-reduction implementation. The npm package metadata declares
MIT but the LICENSE file shipped in the tarball is **Apache-2.0**, so we
comply with the more conservative Apache-2.0 terms (preserve copyright /
attribution; mark modifications, of which there are none).

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

The full Apache-2.0 license text is reproduced in `THIRD_PARTY_LICENSES.txt`,
available at the URL above and in `node_modules/umap-js/LICENSE`.

### jStat (1.9.x)

Statistical-distribution library used for ANOVA, hypothesis testing, and
probabilistic scoring.

- Upstream: https://github.com/jstat/jstat
- License: **MIT**

```
Copyright (c) 2013 jStat

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

### math.js (mathjs, 15.x)

Extensive math library — matrix operations, expression parsing, complex
numbers — used inside the linear-method/PLS pipelines.

- Upstream: https://mathjs.org/ — https://github.com/josdejong/mathjs
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

The full Apache-2.0 license text is reproduced in `THIRD_PARTY_LICENSES.txt`
and available at the URL above. Upstream `node_modules/mathjs/LICENSE` does not
ship a NOTICE file. Copyright is held by *Jos de Jong and contributors*.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into EDA's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |
| wu.js     | 2.x     | MIT        | https://github.com/fitzgen/wu.js        |

`@webgpu/types` (BSD-3-Clause) is a TypeScript types-only dev dependency and
is not present in the runtime bundle.

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published EDA plugin and impose no obligation on users of the
plugin.

`source-map-loader`, `worker-loader`, `css-loader`, `style-loader`,
`ts-loader`, `file-loader` (all MIT) and `@webgpu/types` (BSD-3-Clause) are
build-time only and are not included in the runtime artifact.
