# SureChEMBL — Third-Party Libraries

The `@datagrok/surechembl` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

This package is a thin client for a locally-deployed SureCHEMBL Postgres
database (with the RDKit cartridge enabled). It contains no bundled
third-party JavaScript runtime libraries beyond what the Datagrok platform
already provides as webpack externals; the cheminformatics work is delegated
to the platform-managed `@datagrok/chem` plugin and to the Postgres database
in the Docker container described below.

---

## 1. Bundled in the published artifact (`dist/`)

None. All runtime dependencies are either listed in Section 2 (provided by the
platform) or are Datagrok-internal libraries (`@datagrok-libraries/utils`,
`@datagrok-libraries/db-explorer`) covered by the repo-wide `LICENSE.md`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into SureChEMBL's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component                        | Version | License      | Upstream                                          |
|----------------------------------|---------|--------------|---------------------------------------------------|
| OpenChemLib JS                   | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js        |
| RxJS                             | 6.x     | Apache-2.0   | https://github.com/ReactiveX/rxjs                 |
| cash-dom                         | 8.x     | MIT          | https://github.com/fabiospampinato/cash           |
| Day.js                           | 1.x     | MIT          | https://github.com/iamkun/dayjs                   |
| wu.js                            | 2.x     | MIT          | https://github.com/fitzgen/wu.js                  |

---

## 3. Fetched at runtime from third-party services

None.

---

## 4. Docker container images (`dockerfiles/`)

The SureChEMBL package ships a single Docker image built from
`datagrok/demo_db_surechembl:2024.1`, which packages the public SureCHEMBL
patent-chemistry database into a Postgres / RDKit-cartridge image:

| Component  | Source                                              | License                                                                |
|------------|-----------------------------------------------------|------------------------------------------------------------------------|
| SureChEMBL | https://www.surechembl.org/                          | CC BY-SA 3.0 (database content, EBI)                                  |
| PostgreSQL | https://www.postgresql.org/                          | PostgreSQL License (BSD-style)                                        |
| RDKit      | https://www.rdkit.org/ (Postgres cartridge `rdkit`) | BSD-3-Clause                                                          |

The SureCHEMBL chemistry data is provided by EMBL-EBI and is licensed
**CC BY-SA 3.0**; redistribution of any derived database dump must preserve
attribution to EMBL-EBI and continue to be made available under the same
license.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published SureChEMBL plugin and impose no obligation on users
of the plugin.

The dev dependency on `@datagrok/chem` is MIT — covered by the repo-wide
`LICENSE.md` and detailed in [`../Chem/CREDITS.md`](../Chem/CREDITS.md).
