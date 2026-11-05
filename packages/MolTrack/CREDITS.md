# Mol Track — Third-Party Libraries

The `@datagrok/moltrack` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All components bundled into the published artifact are under permissive
licenses. No copyleft (GPL/LGPL/MPL) component is bundled into the published
Mol Track plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### Crux WebAssembly modules (`vendor/crux/`)

`crux_wasm` (structure standardization and search) and `inchi_wasm` (InChI),
built from [crux-core](https://github.com/datagrok-ai/crux-core) and ported from
RDKit (BSD-3-Clause). Their notices are in
[`vendor/crux/NOTICE.crux-wasm`](vendor/crux/NOTICE.crux-wasm) and
[`vendor/crux/NOTICE.inchi-wasm`](vendor/crux/NOTICE.inchi-wasm).

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Mol Track's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |

---

## 3. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Mol Track plugin and impose no obligation on users of
the plugin.
