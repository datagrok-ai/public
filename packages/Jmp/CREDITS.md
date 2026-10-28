# JMP Reader — Third-Party Libraries

The `@datagrok/jmp-reader` package is distributed under the MIT license that
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

### jmp-convert

The column-property layout in `src/formats/jmp/jmp-properties.ts` (the
property list after a column's descriptor bytes, and the value-label
payload) follows jmp-convert's reader. No source is copied; the layout is
re-implemented in TypeScript.

- Upstream: https://github.com/davidmathers/jmp-convert
- Read at commit: `5024c5aba8be0cbfaa21e5a3be85cd9f8a33c3fc` (2026-09-25)
- License: **MIT**, Copyright (c) 2026 David Mathers

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

---

## 5. Test data (`files/public/`, not bundled)

Data files that researchers published openly, renamed after their source. Used only by the
package tests. CC0 asks for no attribution; those files are credited here anyway. The CC BY 4.0
file is used under that licence (https://creativecommons.org/licenses/by/4.0/), unchanged apart
from its name.

| File | Original name | Licence | Dataset |
|---|---|---|---|
| `dryad-t76hdr7xq-hadmsc.jmp` | `hAdMSC_siTcad+APN.jmp` | CC0 1.0 | Nakamura, Kita, Tanaka et al. (2020), *Adiponectin stimulates exosome release to enhance mesenchymal stem cell driven therapy of heart failure in mice*. Dryad, https://doi.org/10.5061/dryad.t76hdr7xq |
| `dryad-t76hdr7xq-tac2wks.jmp` | `TAC2wks_hMSC(siTcad).jmp` | CC0 1.0 | same dataset |
| `dryad-60nn7-keiser-2016.jmp` | `Keiser et al 2016 - Data File for Proceedings B.jmp` | CC0 1.0 | Keiser, Pinter-Wollman, Augustine et al. (2016), *Data from: Individual differences in boldness influence patterns of social interactions and the transmission of cuticular bacteria among group-mates*. Dryad, https://doi.org/10.5061/dryad.60nn7 |
| `dryad-6b340-acoustic.jmp` | `Acoustic_Dryad.jmp` | CC0 1.0 | Warwick, Travis, Lemmon (2015), *Data from: Geographic variation in the Pine Barrens Treefrog (Hyla andersonii): concordance of genetic, morphometric, and acoustic signal data*. Dryad, https://doi.org/10.5061/dryad.6b340 |
| `dataverse-jylsxj-rlrva.jmp` | `RLRVAsummary_with_code.jmp` | CC0 1.0 | Altman, McDonald (2013), *Replication data for: A Half-Century of Virginia Redistricting Battles: Shifting from Rural Malapportionment to Voting Rights to Public Participation*. Harvard Dataverse, https://doi.org/10.7910/DVN/JYLSXJ |
| `dataverse-hqaihu-figure4.jmp` | `Figure4Replication.jmp` | CC0 1.0 | Altman, McDonald (2015), *Replication Data for: Redistricting and Polarization*. Harvard Dataverse, https://doi.org/10.7910/DVN/HQAIHU |
| `dataverse-jvinzu-rat.jmp` | `Rat.jmp` | CC0 1.0 | Bender (2025), *Closed-Loop Electrical Block of Vagus Nerve Scales from Rodent to Porcine Cardiac Models*. Harvard Dataverse, https://doi.org/10.7910/DVN/JVINZU |
| `figshare-15044331-das28.jmp` | `DAS280M v BL IGS.jmp` | CC BY 4.0 | Thirkettle (2021), *On A Knife Edge: the role of interferon pathways in the adaptive immune response and potential progression to RA pathogenesis*. Newcastle University, figshare, https://doi.org/10.25405/data.ncl.15044331.v1 |
