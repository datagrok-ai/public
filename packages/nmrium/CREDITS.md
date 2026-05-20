# NMRium — Third-Party Libraries

The `@datagrok/nmrium` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

This plugin is a thin wrapper around the upstream **NMRium** React component
and its cheminfo-org / Zakodium ecosystem. The bulk of the bundled code is
third-party and is enumerated below.

All runtime JavaScript dependencies bundled into the published artifact are
under permissive licenses (MIT, Apache-2.0, BSD-3-Clause, CC-BY-4.0). No
copyleft (GPL/LGPL/MPL) JavaScript component is bundled into the published
NMRium plugin.

---

## 1. Bundled in the published artifact (`dist/`)

### NMRium (0.59.x)

The full NMR data viewer / processor React component — the primary feature
of this plugin.

- Upstream: https://www.nmrium.org — https://github.com/cheminfo/nmrium
- License: **MIT**

```
MIT License

Copyright (c) 2021 ChemInfo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### nmr-load-save (0.37.x), nmr-processing (12.x), filelist-utils (1.x), fifo-logger (1.x), openchemlib-utils (6.x)

cheminfo-org NMR / file / chemistry helper libraries used transitively (and
some directly) by NMRium.

- Upstream: https://github.com/cheminfo
- License: **MIT** (each carries an analogous notice)

```
MIT License

Copyright (c) cheminfo

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### OpenChemLib JS (8.x)

Cheminformatics rendering and structure handling. (Although the platform
serves a different OpenChemLib build via webpack externals, NMRium is built
against its own pinned `openchemlib` 8.x and bundles parts of it.)

- Upstream: https://github.com/cheminfo/openchemlib-js
- License: **BSD-3-Clause**

```
Copyright (c) 2015-2017, cheminfo

All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### react-science (6.x)

Zakodium's React-component library used by NMRium for UI scaffolding.

- Upstream: https://github.com/zakodium-oss/react-science
- License: **MIT**

```
MIT License

Copyright (c) 2021 Zakodium

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### react-icons (5.x)

Icon-pack aggregator. Code is **MIT** licensed; individual icon glyphs are
re-distributed under their original licenses (mostly MIT or CC-BY).

- Upstream: https://github.com/react-icons/react-icons
- Code license: **MIT**

```
Copyright 2018 kamijin_fanta <kamijin@live.jp>

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
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
```

Icon-pack licenses (a partial list, as advertised by the upstream
`react-icons/LICENSE`):

- Font Awesome 5/6 — CC BY 4.0
- Material Design icons — Apache-2.0
- Ionicons, Octicons, Feather, Heroicons, Tabler, Bootstrap Icons, … — MIT
- Lucide — ISC
- Simple Icons — CC0 1.0
- Circum Icons — MPL-2.0 (icon glyphs only — no source-code copyleft impact)
- Weather Icons — SIL OFL 1.1
- VS Code Icons / IcoMoon Free — CC BY 4.0
- Typicons / Game Icons — CC BY-SA 3.0 / CC BY 3.0
- Remix Icon / Grommet-Icons — Apache-2.0

Only icons that NMRium actually imports end up in `dist/`; the pack as a whole
is tree-shaken at bundle time.

### @blueprintjs/core (5.x), @blueprintjs/icons (5.x)

Palantir's Blueprint UI toolkit, used by NMRium.

- Upstream: https://blueprintjs.com/
- Copyright: Palantir Technologies, Inc.
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy
> of the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

The full Apache-2.0 license text ships at
`node_modules/@blueprintjs/core/LICENSE`.

### @emotion/react (11.x)

CSS-in-JS library used by NMRium / Blueprint.

- Upstream: https://emotion.sh/
- License: **MIT**

```
MIT License

Copyright (c) Emotion team and other contributors

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
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
```

### cheminfo-font (1.13.x)

Custom icon font assembled by cheminfo from many upstream icon packs. Icon
glyphs are re-distributed under **CC-BY-4.0**.

- Upstream: https://github.com/cheminfo/font
- License: **CC-BY-4.0**

> Licensed under the Creative Commons Attribution 4.0 International License.
> See https://creativecommons.org/licenses/by/4.0/ for the full license text.

The full CC-BY-4.0 text ships at `node_modules/cheminfo-font/LICENSE`.

### React (18.x), React-DOM (18.x)

Required peer of NMRium and react-science. Bundled as part of the NMRium plugin.

- Upstream: https://reactjs.org/
- Copyright: Meta Platforms, Inc. and affiliates
- License: **MIT**

```
MIT License

Copyright (c) Meta Platforms, Inc. and affiliates.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Other transitive bundled libraries

NMRium pulls in dozens of additional small permissive dependencies (e.g.
ml-* packages, d3-*, lodash, plotly subset components, jsondiffpatch). All
were verified at integration time to be MIT, BSD, ISC, or Apache-2.0. Their
individual notices are preserved verbatim by webpack's
`license-webpack-plugin` output and by the unmodified runtime modules in
`dist/`.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into NMRium's `dist/` — they are provided
once by the platform host and shared across all packages.

| Component      | Version | License      | Upstream                                |
|----------------|---------|--------------|-----------------------------------------|
| OpenChemLib JS | 7.x     | BSD-3-Clause | https://github.com/cheminfo/openchemlib-js |
| cash-dom       | 8.x     | MIT          | https://github.com/fabiospampinato/cash |
| Day.js         | 1.x     | MIT          | https://github.com/iamkun/dayjs         |

OpenChemLib JS is also referenced by NMRium's `package.json`. The upstream
`openchemlib-js` npm package is distributed under BSD-3-Clause:
*Copyright (c) 2015-2017, cheminfo.*

---

## 3. Docker container images (`dockerfiles/`)

None. NMRium does not ship a Docker image.

---

## 4. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published NMRium plugin and impose no obligation on users of
the plugin.
