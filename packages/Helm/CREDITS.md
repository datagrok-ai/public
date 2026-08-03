# Helm — Third-Party Libraries

The `@datagrok/helm` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, ISC, BSD-3-Clause / Academic Free License v2.1).
No copyleft (GPL/LGPL/MPL) component is bundled into the published Helm
plugin artifact.

The Pistoia HELM Web Editor and Scilligence JSDraw2 sources, which were
historically vendored under `helm/JSDraw/`, are now consumed via
`@datagrok-libraries/helm-web-editor` (an internal Datagrok library covered by
the repo-wide MIT umbrella). Their original upstream licenses are reproduced
below for completeness because the vendored copies remain in the repository.

---

## 1. Bundled in the published artifact (`dist/`)

### Dojo Toolkit (1.10.10)

The full Dojo Toolkit (`dojo`, `dijit`, `dojox`) is vendored under
`vendor/dojo-1.10.10/` and webpack-bundled into a separate `package-dojo.js`
chunk via `require.context()`. Required by JSDraw2 / HELM Web Editor.

- Upstream: https://dojotoolkit.org/ — https://github.com/dojo/dojo
- License: **BSD-3-Clause** OR **Academic Free License v2.1** (dual-licensed,
  at the user's choice — Datagrok elects the BSD-3-Clause terms)

```
Dojo is available under *either* the terms of the BSD 3-Clause "New" License
*or* the Academic Free License version 2.1. As a recipient of Dojo, you may
choose which license to receive this code under.

Copyright (c) 2005-2015, The Dojo Foundation
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
  * Neither the name of the Dojo Foundation nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

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

**Modification note.** `core/client/d4`-style monkey patches are applied to
Dojo at runtime (not in the vendored source) — see
`HelmPackage.initHELMWebEditor()` in `src/package-utils.ts`, which patches
`dojox.gfx.svg.Text.prototype.getTextWidth` to prevent an infinite loop. The
vendored Dojo files themselves are byte-identical to the upstream 1.10.10
release.

### Scilligence JSDraw2 Lite & Pistoia HELM Web Editor (vendored copy under `helm/JSDraw/`)

Historically vendored uncompressed sources of the JSDraw2 Lite molecule
sketcher (Scilligence) and the Pistoia HELM Web Editor. The active code path
loads `@datagrok-libraries/helm-web-editor` (the internal Datagrok library)
rather than these files; they are kept in the repository as a reference and
fallback. If a build configuration ever bundles them again, the following
notice applies.

- Upstream (JSDraw2 Lite): https://github.com/scilligence/JSDraw.Lite
- Upstream (HELM Web Editor): https://github.com/PistoiaHELM/HELMWebEditor
- License: **MIT** (both projects)

```
The MIT License (MIT)

Copyright (c) Scilligence Corporation (JSDraw2 Lite) and the Pistoia Alliance
(HELM Web Editor)

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

### lru-cache (10.4.x)

LRU cache used to cap the rendered-monomer / off-screen-editor caches in the
HELM cell renderer.

- Upstream: https://github.com/isaacs/node-lru-cache
- License: **ISC**

```
The ISC License

Copyright (c) 2010-2023 Isaac Z. Schlueter and Contributors

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
```

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Helm's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component       | Version | License      | Upstream                                   |
|-----------------|---------|--------------|--------------------------------------------|
| OpenChemLib JS  | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js |
| RxJS            | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs          |
| cash-dom        | 8.x     | MIT          | https://github.com/fabiospampinato/cash    |
| Day.js          | 1.x     | MIT          | https://github.com/iamkun/dayjs            |
| wu.js           | 2.x     | MIT          | https://github.com/fitzgen/wu.js           |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Helm plugin and impose no obligation on users of the
plugin.

The peer/devDependencies on other Datagrok plugins (`@datagrok/bio`,
`@datagrok/chem`) are MIT — covered by the repo-wide `LICENSE.md`. Each of
those plugins maintains its own `CREDITS.md` for its own third-party content.
