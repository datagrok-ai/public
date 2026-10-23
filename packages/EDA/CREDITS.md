# EDA — Third-Party Libraries

The `@datagrok/eda` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, Apache-2.0, BSD-3-Clause, BSL-1.0, Unlicense), except
the Eigen port inside the Crux module, which is MPL-2.0 (file-level copyleft,
see below). No GPL/LGPL component is bundled into the published EDA plugin
artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### Crux (`crux/eda/crux_wasm_bg.wasm`)

Datagrok's Crux engine (crux-core, the commit pinned in `crux/engine.toml`),
built as the `crux-wasm` crate with the feature list of `crux/engine.toml` and
shipped as `crux/eda/crux_wasm.js` + `crux_wasm_bg.wasm` (copied to
`dist/crux/eda/`). It runs PCA, PLS, linear and logistic regression, XGBoost
and SVM. The module contains Rust ports of the following upstreams (the crates
compiled into the shipped binary; licenses as declared by the crux crates:
`Cargo.toml` `license`, the `LICENSE*` files and `LICENSE` consts of the
crates, the model `upstreams()`, the pins in `upstreams/*.toml`):

| Upstream | Version | License | Crux crate |
|---|---|---|---|
| scikit-learn (incl. its vendored, modified LIBSVM 3.10 and LIBLINEAR) | 1.9.1 | BSD-3-Clause | `ml` |
| astartes (QSAR splits) | 1.3.3 | MIT | `ml` |
| SIMPD (QSAR splits) | 403f758 | MIT | `ml` |
| pymoo (SIMPD's optimizer) | 0.5.0 | Apache-2.0 | `ml` |
| XGBoost | 3.4.1 | Apache-2.0 | `xgboost` |
| LightGBM (compiled in as a dependency of `shap`; no LightGBM functions in the module) | 4.7.0 | MIT | `lightgbm` |
| CatBoost | 1.2.10 | Apache-2.0 | `catboost` |
| imbalanced-learn | 0.14.2 | MIT | `imblearn` |
| shap | 0.53.0 | MIT | `shap` |
| crepes | 0.9.1 | BSD-3-Clause | `conformal` |
| venn-abers | 1.5.4 | MIT | `conformal` |
| MAPIE (subset) | 1.5.0 | BSD-3-Clause | `conformal` |
| Optuna | 5.0.0 | MIT | `hpo` |
| cmaes | 0.13.1 | MIT | `hpo` |
| bayesian-optimization (bayes_opt) | 3.4.0 | MIT | `hpo` |
| SciPy | 1.18.1 | BSD-3-Clause | `scipy` |
| Boost.Math (as vendored by SciPy 1.18.1) | c114cf4 | BSL-1.0 | `scipy` |
| Reference LAPACK (+ reference BLAS) | 3.12.0 | BSD-3-Clause | `lapack` |
| HiGHS | 1.12.0 | MIT | `highs` |
| NumPy (random / linalg / sort building blocks) | 1.26.4 | BSD-3-Clause | `pycompat` |
| David M. Gay's dtoa.c (float formatting) | | Lucent notice (permissive) | `pycompat` |
| RDKit (molecule model, SMILES / SMARTS, descriptors, fingerprints, MolStandardize, reactions) | 2026.03.6 | BSD-3-Clause | `rk`, `smi-compiler`, `smarts-parser`, `molstandardize`, `rxn` |
| Avalon Toolkit (Avalon fingerprints, as RDKit builds it) | 2.0.5-pre.3 | BSD-3-Clause | `rk` |
| Eigen (dense decompositions, self-adjoint eigensolver, SVD) | 3.4.0 / 3.3.9 | **MPL-2.0** | `eigen` |

- scikit-learn: https://scikit-learn.org — `Copyright (c) 2007-2026 The scikit-learn developers`
- astartes: https://github.com/JacksonBurns/astartes — `Copyright (c) 2023 Jackson Burns`
- SIMPD: https://github.com/rinikerlab/molecular_time_series — `Copyright (c) 2022 Greg Landrum`
- pymoo: https://github.com/anyoptimization/pymoo
- XGBoost: https://github.com/dmlc/xgboost — `Copyright (c) 2019 by Contributors`
- LightGBM: https://github.com/microsoft/LightGBM
- CatBoost: https://github.com/catboost/catboost
- LIBSVM (inside scikit-learn): https://www.csie.ntu.edu.tw/~cjlin/libsvm/
- imbalanced-learn: https://github.com/scikit-learn-contrib/imbalanced-learn — `Copyright (c) 2014-2020 The imbalanced-learn developers`
- shap: https://github.com/shap/shap
- crepes: https://github.com/henrikbostrom/crepes — `Copyright 2026 Henrik Boström`
- venn-abers: https://github.com/ip200/venn-abers — `Copyright (c) 2023 Ivan Petej`
- MAPIE: https://github.com/scikit-learn-contrib/MAPIE — `Copyright (c) 2021, Vianney TAQUET, Grégoire MARTINON, Nicolas BRUNEL, Issam IBNOUHSEIN, François DEHEEGER and MAPIE contributors`
- Optuna: https://github.com/optuna/optuna — `Copyright (c) 2018 Preferred Networks, Inc.`
- cmaes: https://github.com/CyberAgentAILab/cmaes — `Copyright (c) 2020 CyberAgent, Inc.`
- bayesian-optimization: https://github.com/bayesian-optimization/BayesianOptimization — `Copyright (c) 2014 Fernando M. F. Nogueira`
- SciPy: https://github.com/scipy/scipy — `Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers`
- Boost.Math: https://github.com/boostorg/math
- LAPACK: https://github.com/Reference-LAPACK/lapack — `Copyright (c) 1992-2023 The University of Tennessee and The University of Tennessee Research Foundation; Copyright (c) 2000-2023 The University of California Berkeley; ...`
- HiGHS: https://github.com/ERGO-Code/HiGHS
- NumPy: https://github.com/numpy/numpy — `Copyright (c) 2005-2023, NumPy Developers`
- dtoa.c: `Copyright (c) 1991, 2000, 2001 by Lucent Technologies` (author David M. Gay)
- RDKit: https://github.com/rdkit/rdkit — `Copyright (C) Greg Landrum and other RDKit contributors`
- Avalon Toolkit: https://github.com/rdkit/ava-formake — `Copyright 2001-2011 Novartis Pharma AG`
- Eigen: https://gitlab.com/libeigen/eigen — MPL-2.0 (https://mozilla.org/MPL/2.0/). The
  Source Code Form of the ported Eigen files is the `eigen` crate of crux-core at the commit
  pinned in `crux/engine.toml`; the original is available at the URL above.

The Apache-2.0 license text is reproduced in `THIRD_PARTY_LICENSES.txt`. The
BSD-3-Clause notice of LIBSVM, which scikit-learn's port carries:

```
Copyright (c) 2000-2009 Chih-Chung Chang and Chih-Jen Lin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither name of copyright holders nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### Apache Arrow JavaScript (apache-arrow, 21.x)

Arrow IPC encoding of the data exchanged with the Crux module.

- Upstream: https://github.com/apache/arrow-js
- License: **Apache-2.0**; NOTICE: `Apache Arrow JavaScript, Copyright
  2017-2025 The Apache Software Foundation. This product includes software
  developed at The Apache Software Foundation (http://www.apache.org/).`

### fast-sha256 (1.3.x)

SHA-256 of a model blob (the Crux model cache key).

- Upstream: https://github.com/dchest/fast-sha256-js
- License: **Unlicense** (public domain)

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
