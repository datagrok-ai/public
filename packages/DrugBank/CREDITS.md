# DrugBank — Third-Party Libraries

The `@datagrok/drug-bank` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

DrugBank does not bundle any third-party JavaScript runtime libraries beyond the
platform-provided externals listed below. All chemistry operations are delegated
at runtime to the `@datagrok/chem` plugin; no copyleft (GPL/LGPL/MPL) component
is bundled into the published artifact.

The bundled drug dataset (`files/drugbank-open-structures.d42`) is derived from
the [DrugBank Open Structures](https://go.drugbank.com/releases/latest#open-data)
release, which is distributed under the Creative Commons Attribution-NonCommercial
4.0 International License (CC BY-NC 4.0). Attribution: © OMx Personal Health
Analytics, Inc. — https://go.drugbank.com/. Users redistributing this dataset
outside the Datagrok platform must comply with the upstream CC BY-NC 4.0 terms.

---

## 1. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into DrugBank's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component       | Version | License      | Upstream                                   |
|-----------------|---------|--------------|--------------------------------------------|
| OpenChemLib JS  | 8.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js |
| cash-dom        | 8.x     | MIT          | https://github.com/fabiospampinato/cash    |
| Day.js          | 1.x     | MIT          | https://github.com/iamkun/dayjs            |

OpenChemLib JS is also listed in DrugBank's `package.json` `sources`. The
upstream `openchemlib-js` npm package is distributed under BSD-3-Clause:
*Copyright (c) 2015-2017, cheminfo.*

---

## 2. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published DrugBank plugin and impose no obligation on users of
the plugin.
