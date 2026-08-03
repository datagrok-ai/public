# Boltz1 — Third-Party Libraries

The `@datagrok/boltz1` package is distributed under the MIT license that
covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT). No copyleft (GPL/LGPL/MPL) component is
bundled into the published Boltz1 plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### js-yaml (4.x)

YAML 1.2 parser/dumper used for reading/writing Boltz-1 input configuration.

- Upstream: https://github.com/nodeca/js-yaml
- License: **MIT**

```
(The MIT License)

Copyright (C) 2011-2015 by Vitaly Puzrin

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

These libraries are not bundled into Boltz1's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |
| cash-dom  | 8.x     | MIT        | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

---

## 3. Docker container (`dockerfiles/boltz`) — Boltz-1 inference service

The Boltz1 package ships a Docker image (`dockerfiles/boltz/Dockerfile`)
that runs the upstream **Boltz-1** structure-prediction model as a Flask
service. The image is built from `nvidia/cuda:11.8.0-runtime-ubuntu22.04`
and pulls third-party software at build time via micromamba and a
git-clone of the upstream Boltz repository:

| Component                            | Source                                | License                                |
|--------------------------------------|---------------------------------------|----------------------------------------|
| Boltz-1                              | https://github.com/jwohlwend/boltz    | MIT                                    |
| NVIDIA CUDA runtime, cuDNN           | nvidia/cuda Docker base image         | NVIDIA EULA — redistribution restricted |
| micromamba                           | https://micro.mamba.pm/               | BSD-3-Clause                           |
| Python 3.10                          | https://www.python.org/               | PSF License                            |
| Flask, flask-cors                    | conda-forge                           | BSD-3-Clause                           |
| PyTorch, torchvision, torchaudio     | conda pytorch channel                 | BSD-3-Clause-style (PyTorch license)   |
| AWS CLI                              | conda-forge (build-time only)         | Apache-2.0                             |

The image also pulls pre-trained Boltz-1 model weights from the public
`s3://datagrok-data/models/boltz` bucket at build time.

If the Boltz1 image is published to a public registry, the NVIDIA CUDA EULA
should be re-reviewed for the redistribution scenario.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published Boltz1 plugin and impose no
obligation on users of the plugin.
