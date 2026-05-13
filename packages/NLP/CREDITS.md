# NLP — Third-Party Libraries

The `@datagrok/nlp` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT, Apache-2.0, BSD-3-Clause). No copyleft
(GPL/LGPL/MPL) JavaScript component is bundled into the published NLP plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### @xenova/transformers (2.17.x)

JavaScript port of the Hugging Face `transformers` library, used to run
text-embedding and similarity models in the browser. Bundles `onnxruntime-web`
under the hood; the WebAssembly artifacts in `dist/*.wasm` are the
ONNX Runtime WASM modules emitted by transformers.js / onnxruntime-web.

- Upstream: https://github.com/xenova/transformers.js
- Author: Joshua (Xenova)
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy
> of the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
> License for the specific language governing permissions and limitations
> under the License.

The full Apache-2.0 license text ships at `node_modules/@xenova/transformers/LICENSE`.

### onnxruntime-web (1.14.x) — transitive via @xenova/transformers

ONNX Runtime compiled to WebAssembly. The two `.wasm` files emitted into
`dist/` (`8473fcbfb6e85ca6c852.wasm`, `9a8fbf37666e32487835.wasm`) are
ONNX Runtime kernels.

- Upstream: https://github.com/microsoft/onnxruntime
- Copyright: Microsoft Corporation
- License: **MIT**

```
MIT License

Copyright (c) Microsoft Corporation

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

### aws-sdk (2.x)

AWS SDK for JavaScript, used to integrate with Amazon services (e.g. S3,
Comprehend) for cloud-based NLP workflows.

- Upstream: https://github.com/aws/aws-sdk-js
- License: **Apache-2.0**

```
AWS SDK for JavaScript
Copyright 2012-2017 Amazon.com, Inc. or its affiliates. All Rights Reserved.

This product includes software developed at
Amazon Web Services, Inc. (http://aws.amazon.com/).
```

> Licensed under the Apache License, Version 2.0. The full license text ships
> at `node_modules/aws-sdk/LICENSE.txt` and is available at
> http://www.apache.org/licenses/LICENSE-2.0.

### umap-js (1.4.x)

UMAP dimensionality reduction used for embedding-space visualization. The
npm package metadata declares MIT but the LICENSE file shipped in the tarball
is **Apache-2.0**, so we comply with the more conservative Apache-2.0 terms
(preserve copyright/attribution; mark modifications, of which there are none).

- Upstream: https://github.com/PAIR-code/umap-js
- Author: Andy Coenen (Google PAIR)
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy
> of the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

### stemmer (2.0.x)

Porter stemming algorithm — used for keyword extraction.

- Upstream: https://words.github.io/stemmer/
- Author: Titus Wormer
- License: **MIT**

```
(The MIT License)

Copyright (c) 2014 Titus Wormer <tituswormer@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
'Software'), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

### util (0.12.x)

Browser-compatible polyfill for Node's `util` module.

- Upstream: https://github.com/browserify/node-util
- License: **MIT**

```
Copyright Joyent, Inc. and other Node contributors. All rights reserved.
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
```

### @webgpu/types

WebGPU TypeScript type definitions. Type-only — no runtime code is bundled.

- Upstream: https://github.com/gpuweb/types
- License: **BSD-3-Clause**

```
Copyright 2022 WebGPU Developers

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into NLP's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License    | Upstream                          |
|-----------|---------|------------|-----------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs |

---

## 3. Fetched at runtime from third-party services (not bundled)

The `@xenova/transformers` runtime fetches model weights on first use from
the **Hugging Face Hub** (https://huggingface.co/) by default. Individual
models have their own licenses (commonly Apache-2.0 or MIT for sentence
transformers). No model weights are bundled with this plugin.

---

## 4. Docker container images (`dockerfiles/`)

### `dockerfiles/` — Sentence-Transformers embedding service

The image is built `FROM python:3.11-slim` and pre-downloads the
`all-MiniLM-L6-v2` sentence-transformer model at build time so the service
starts offline.

| Component                  | Source                                                             | License      |
|----------------------------|--------------------------------------------------------------------|--------------|
| Flask                      | https://flask.palletsprojects.com/                                 | BSD-3-Clause |
| sentence-transformers      | https://github.com/UKPLab/sentence-transformers                    | Apache-2.0   |
| transformers (Hugging Face)| https://github.com/huggingface/transformers (transitive)           | Apache-2.0   |
| PyTorch                    | https://pytorch.org/ (transitive)                                  | BSD-3-style  |
| NumPy                      | https://numpy.org/                                                 | BSD-3-Clause |
| `all-MiniLM-L6-v2` weights | https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2      | Apache-2.0   |
| Python stdlib              | https://www.python.org/                                            | PSF          |
| Debian slim base packages  | https://www.debian.org/                                            | various      |

The Debian base image includes GNU/Linux system libraries under copyleft
licenses (glibc — LGPL; coreutils, bash, apt — GPL). These are standard OS
components used as-is and are subject to the obligations of their upstream
licenses if the image itself is redistributed.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published NLP plugin and impose no obligation on users of the
plugin.
