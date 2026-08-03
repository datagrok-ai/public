# Metabolic Graph — Third-Party Libraries

The `@datagrok/metabolicgraph` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

**Important licensing note:** This package bundles `glpk.js`, which is licensed
under **GPL-3.0** (a strong copyleft license). As a consequence, the published
Metabolic Graph plugin artifact is effectively redistributable under GPL-3.0
terms — copyleft obligations apply to anyone who redistributes the bundled
plugin. The rest of the Metabolic Graph source code remains MIT-licensed under
the repo-wide `LICENSE.md`.

---

## 1. Bundled in the published artifact (`dist/`)

### Escher (1.7.3)

Web application for building, sharing, and embedding data-rich visualizations
of metabolic pathways. A pre-built copy of `escher.js` is checked in at
`escher_dist/escher.js`; the upstream source is mirrored under `escher_src/`.

- Upstream: https://escher.github.io — https://github.com/zakandrewking/escher
- Author: Zachary King
- License: **MIT**

```
The MIT License (MIT)

This software is Copyright © 2015 The Regents of the University of
California. All Rights Reserved.

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

### glpk.js (4.0.2)

Emscripten port of the GNU Linear Programming Kit (GLPK), used for
flux-balance analysis (FBA) and related linear-programming problems on
metabolic models. The compiled `glpk.all.wasm` is bundled into `dist/`.

- Upstream: https://github.com/jvail/glpk.js
- Author: Jan Vaillant (port); GLPK by the GNU Project
- License: **GPL-3.0** (copyleft — see top of file)

> Licensed under the GNU General Public License, Version 3. See
> https://www.gnu.org/licenses/gpl-3.0.html for the full license text.
> The full license text also ships at `node_modules/glpk.js/LICENSE`.
>
> You should have received a copy of the GNU General Public License along
> with this program. If not, see http://www.gnu.org/licenses/.

### COBRA / Eigen-based sampler (`src/cobra/sampler.wasm`)

Custom sampler for COBRA-style flux sampling, written in C++ and compiled to
WebAssembly via Emscripten. Source: `src/cobra/sampler-minimal.cpp`. The
compiled `sampler.wasm` is bundled into `dist/`.

- Original code: Datagrok contributors — covered by the repo-wide MIT license
- Linked header library: **Eigen 3.4.0** (linear algebra, header-only)

Eigen is primarily licensed under the **MPL-2.0** (Mozilla Public License, v2),
with some modules under BSD or LGPL. The MPL-2.0 is a weak file-level copyleft
that does not impose obligations on the rest of the bundle.

- Upstream: https://eigen.tuxfamily.org/
- License: **MPL-2.0** (with some BSD/LGPL files)

> Mozilla Public License, version 2.0. See https://www.mozilla.org/MPL/2.0/
> for the full license text.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Metabolic Graph's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Docker container images (`dockerfiles/`)

None. Metabolic Graph does not ship a Docker image.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Metabolic Graph plugin and impose no obligation on
users of the plugin.
