# PreclinicalCase — Third-Party Libraries

The `preclinicalcase` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

The TypeScript bundle has only one third-party runtime JavaScript dependency
(`x2js`, Apache-2.0). The bundled CDISC Open Rules Engine (CORE) binary used
for SEND validation is fetched at container-build time and runs server-side
inside the Docker image (Section 4).

---

## 1. Bundled in the published artifact (`dist/`)

### x2js (3.4.x)

XML ↔ JavaScript object converter used to parse SEND study configuration
XML files.

- Upstream: https://github.com/abdmob/x2js
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

The full Apache-2.0 license text is available at the URL above and in
`node_modules/x2js/LICENSE`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into PreclinicalCase's `dist/` — they are
provided once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                                |
|-----------|---------|------------|-----------------------------------------|
| Day.js    | 1.x     | MIT        | https://github.com/iamkun/dayjs         |
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs       |

---

## 3. Fetched at runtime from third-party CDNs (not bundled)

None.

---

## 4. Docker container images (`dockerfiles/`)

These images are built and distributed separately from the JavaScript bundle.

### `dockerfiles/` — CDISC SEND validation Celery worker

The image is built from `python:3.12-slim` and packages a Celery worker that
executes SEND-data validation against the CDISC SENDIG ruleset. The
`core-ubuntu-latest` binary is downloaded from the upstream cdisc-rules-engine
GitHub release at image-build time (not at deploy time).

| Component | Source | License |
|-----------|--------|---------|
| Python 3.12 | https://www.python.org/ | PSF License |
| Celery | https://github.com/celery/celery | BSD-3-Clause |
| filelock | https://github.com/tox-dev/py-filelock | Unlicense (public domain) |
| datagrok-celery-task | (Datagrok-internal) | covered by repo MIT |
| datagrok-api (Python) | (Datagrok-internal) | covered by repo MIT |
| **CDISC Open Rules Engine (`core`)** | https://github.com/cdisc-org/cdisc-rules-engine (release `core-ubuntu-latest.zip`) | **MIT** |

The `core` binary is the reference engine maintained by CDISC and the Open
Rules Initiative; it is redistributed inside the container under MIT terms.
The container also installs Debian system libraries (`curl`, `unzip`, `jq`,
`ca-certificates`) provided by the upstream Debian distribution under their
respective licenses.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI and the TypeScript / webpack toolchain. These are
**not** redistributed as part of the published PreclinicalCase plugin and
impose no obligation on users of the plugin.
