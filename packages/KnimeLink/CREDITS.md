# KnimeLink — Third-Party Libraries

The `@datagrok/knimelink` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
All runtime dependencies are provided by the Datagrok platform as webpack
externals — nothing third-party is bundled into the published `dist/`.

---

## 1. Bundled in the published artifact (`dist/`)

None. KnimeLink's own source is the only code in the bundle (apart from the
webpack runtime). All external libraries it uses are platform externals listed
in Section 2.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into KnimeLink's `dist/` — they are provided
once by the platform host and shared across all packages. Listed here for full
transparency; their license texts ship with the platform itself.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party services (not bundled)

KnimeLink communicates with a user-configured **KNIME Business Hub** instance
(default `https://hub.knime.com`) over its REST API. The KNIME Hub server
software and any workflows hosted on it are operated by the user/customer and
licensed by KNIME AG separately; nothing from KNIME is redistributed by this
package.

---

## 4. Docker container images (`dockerfiles/`)

None. KnimeLink does not ship a Docker image.

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published KnimeLink plugin and impose no obligation on users of
the plugin.
