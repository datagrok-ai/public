# JMP — Third-Party Libraries

The `@datagrok/jmp` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
Its preview is built with `@datagrok-libraries/u2`, part of the same repository.
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

The package has no bundled runtime dependencies. Gzip decompression uses the
browser's built-in `DecompressionStream`. No copyleft (GPL/LGPL/MPL) component
is bundled into the published JMP plugin.

---

## 1. Ported source code

### JMPReader.jl

The `.jmp` parser in `src/formats/jmp/` is a TypeScript port of JMPReader.jl.
The test fixtures in `files/jmpreader/` are copied unchanged from that
repository's `test/` directory.

- Upstream: https://github.com/jaakkor2/JMPReader.jl
- Ported from commit: `9fd68aa2651a8dde002ecdbd06af4d8d832b18f2` (v0.2.4, 2026-10-16)
- License: **MIT**

```
Copyright (C) 2024 Jaakko Ruohio

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

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into JMP's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component | Version | License    | Upstream                          |
|-----------|---------|------------|-----------------------------------|
| RxJS      | 6.x     | Apache-2.0 | https://github.com/ReactiveX/rxjs |

---

## 3. Docker container images (`dockerfiles/`)

None. JMP does not ship a Docker image.

---

## 4. Development-only tools

`tools/jmp2json.jl` runs JMPReader.jl (MIT), JSON3.jl (MIT), DataFrames.jl (MIT)
and ColorTypes.jl (MIT) under Julia to generate the golden reference files in
`files/golden/`. Neither Julia nor these libraries are bundled or required at
runtime. `tools/jmp-corpus.ts` runs under `tsx` (MIT), and the behavioural
features in `bdd/` run under `@datagrok-libraries/bdd` and Playwright
(Apache-2.0); none of these are bundled either.

The `datagrok-tools` CLI used during build/test transitively pulls in
`puppeteer-screen-recorder`, which references `@ffmpeg-installer/ffmpeg`
(LGPL-2.1) and `@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not**
redistributed as part of the published JMP plugin and impose no obligation on
users of the plugin.
