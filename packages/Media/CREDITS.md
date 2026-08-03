# Media — Third-Party Libraries

The `@datagrok/media` package is distributed under the MIT license that covers
the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).

Media has no third-party runtime dependencies — its only dependency is the
internal `datagrok-api` (covered by the repo-wide MIT). The package wraps the
browser's built-in `<audio>` and `<video>` elements as Datagrok file viewers
and bundles no external libraries into the published `dist/`.

---

## 1. Bundled in the published artifact (`dist/`)

None.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

None.

---

## 3. Docker container images (`dockerfiles/`)

None. Media does not ship a Docker image.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Media plugin and impose no obligation on users of
the plugin.
