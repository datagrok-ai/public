# NodeJSDemo — Third-Party Libraries

The `nodetodo` (NodeJSDemo) package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

The TypeScript plugin itself bundles no third-party JavaScript runtime
dependencies — its only `package.json` dependency is the internal
`datagrok-api` (covered by the repo-wide MIT). The `dockerfiles/todo-app/`
directory ships a Node.js demo service whose third-party content is itemized
below.

---

## 1. Bundled in the published artifact (`dist/`)

None.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into NodeJSDemo's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Docker container images (`dockerfiles/`)

### `dockerfiles/todo-app` — Node.js TODO demo service

The image is built `FROM node:20-slim` (multi-stage) and runs a small Express
service that demonstrates a Node-side use of the Datagrok API.

| Component                  | Source                                       | License      |
|----------------------------|----------------------------------------------|--------------|
| Node.js                    | https://nodejs.org/                          | MIT          |
| Express                    | https://expressjs.com/                       | MIT          |
| node-postgres (`pg`)       | https://node-postgres.com/                   | MIT          |
| `datagrok-api` (Node build)| https://github.com/datagrok-ai/public        | MIT          |
| Debian slim base packages  | https://www.debian.org/                      | various      |

The Debian base image includes GNU/Linux system libraries under copyleft
licenses (glibc — LGPL; coreutils, bash, apt — GPL). These are standard OS
components used as-is and are subject to the obligations of their upstream
licenses if the image itself is redistributed.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published NodeJSDemo plugin and impose no obligation on users
of the plugin.
