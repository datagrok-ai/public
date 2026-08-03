# KetcherSketcher — Third-Party Libraries

The `@datagrok/ketcher-sketcher` package is distributed under the MIT license
that covers the rest of the `public/` repository (see
[`../../LICENSE.md`](../../LICENSE.md)). It incorporates the open-source
components listed below; this file reproduces the attribution and notices
required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (Apache-2.0, MIT). No copyleft (GPL/LGPL/MPL) component
is bundled into the published KetcherSketcher plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### Ketcher — `ketcher-core`, `ketcher-react`, `ketcher-standalone`, `ketcher-macromolecules` (3.12.x)

EPAM Life Sciences's open-source 2D molecule sketcher. The four npm packages
implement the data model, React UI, standalone (no-server) mode, and the
macromolecules editor respectively. They are all bundled into KetcherSketcher's
`dist/` to provide the "Ketcher" choice in any Datagrok molecule sketcher
dropdown.

The `ketcher-standalone` distribution embeds **Indigo** (`indigo-ketcher`) as a
WebAssembly module (`indigo-ketcher-1.27.0.wasm` and the matching JS glue) for
server-less molfile / SMILES / SMARTS conversion.

- Upstream: https://lifescience.opensource.epam.com/ketcher/index.html —
  https://github.com/epam/ketcher
- Authors: EPAM Life Sciences
- License: **Apache-2.0** (declared by all four `ketcher-*` `package.json`
  files; the upstream repository ships a top-level Apache-2.0 LICENSE)

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
> License for the specific language governing permissions and limitations
> under the License.

The full Apache-2.0 license text is available at the URL above. The npm
tarballs do not ship a separate NOTICE file. Copyright is held by *EPAM
Systems, Inc.* and contributors.

The four CSS class-name and minor visibility workarounds applied in
`css/editor.css` and the `min-width: 0` adjustment in `src/ketcher.tsx` are
modifications applied via additive CSS / runtime style; they do not modify
the upstream Ketcher source files.

### indigo-ketcher (1.27.0)

Indigo cheminformatics toolkit, compiled to WebAssembly and bundled by
`ketcher-standalone`.

- Upstream: https://lifescience.opensource.epam.com/indigo/api/index.html —
  https://github.com/epam/Indigo
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0

The full Apache-2.0 license text is in
`node_modules/indigo-ketcher/LICENSE`. Copyright is held by *EPAM Systems,
Inc.*

### Static SDF templates

`templates/fg.sdf` and `templates/library.sdf` (functional groups and template
library) are routed through `file-loader` and shipped as static assets. They
are derived from the Ketcher template set (Apache-2.0).

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into KetcherSketcher's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

`react` and `react-dom` are required peer dependencies of `ketcher-react` and
are bundled transitively (MIT, *Copyright (c) Meta Platforms, Inc. and
affiliates*).

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published KetcherSketcher plugin and impose no obligation on
users of the plugin.

The peer/devDependency on `@datagrok/chem` is MIT — covered by the repo-wide
`LICENSE.md`. That plugin maintains its own `CREDITS.md` for its third-party
content.
