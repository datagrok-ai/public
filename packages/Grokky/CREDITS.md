# Grokky — Third-Party Libraries

The `@datagrok/grokky` package is distributed under the MIT license that
covers the rest of the `public/` repository (see [`../../LICENSE.md`](../../LICENSE.md)).
It incorporates the open-source components listed below; this file reproduces
the attribution and notices required by their respective licenses.

All runtime dependencies bundled into the published artifact are under
permissive licenses (MIT, Apache-2.0). No copyleft (GPL/LGPL/MPL) component is
bundled into the published Grokky plugin artifact.

---

## 1. Bundled in the published artifact (`dist/`)

### Vercel AI SDK — `ai`, `@ai-sdk/provider`, `@ai-sdk/anthropic`, `@ai-sdk/openai`, `@ai-sdk/amazon-bedrock`

Provider-agnostic LLM client used by the in-browser AI features (chat panels,
search, prompt routing).

- Upstream: https://ai-sdk.dev/ — https://github.com/vercel/ai
- License: **Apache-2.0**

> Copyright 2023 Vercel, Inc.
>
> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
> License for the specific language governing permissions and limitations
> under the License.

The full Apache-2.0 license text is available at the URL above and in
`node_modules/ai/LICENSE`. The upstream tarball does not ship a NOTICE file.

### openai (6.x)

Official OpenAI JavaScript SDK.

- Upstream: https://github.com/openai/openai-node
- License: **Apache-2.0**

> Licensed under the Apache License, Version 2.0 (the "License"); you may not
> use this file except in compliance with the License. You may obtain a copy of
> the License at http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
> WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

Copyright is held by *OpenAI, Inc.*

### @modelcontextprotocol/sdk (1.x)

Anthropic's official Model Context Protocol SDK — used by the in-browser
client to interact with the MCP server (see Section 4 below).

- Upstream: https://modelcontextprotocol.io/ —
  https://github.com/modelcontextprotocol/typescript-sdk
- License: **MIT**

```
MIT License

Copyright (c) 2024 Anthropic, PBC

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

### @types/dom-chromium-ai

TypeScript-only type declarations for the experimental Chromium on-device AI
APIs. Build-time only — no runtime code is bundled. License: **MIT**.

---

## 2. Linked at runtime via the Datagrok platform (webpack externals)

These libraries are not bundled into Grokky's `dist/` — they are provided once
by the platform host and shared across all packages.

| Component | Version | License | Upstream                                |
|-----------|---------|---------|-----------------------------------------|
| cash-dom  | 8.x     | MIT     | https://github.com/fabiospampinato/cash |
| Day.js    | 1.x     | MIT     | https://github.com/iamkun/dayjs         |

---

## 3. Fetched at runtime from third-party LLM providers (not bundled)

Grokky is configured to send chat / completion / tool-use requests to one of
several LLM providers (Anthropic Claude, OpenAI, AWS Bedrock, etc.) selected
by the deployed Datagrok server's configuration. These provider services are
governed by their own respective terms of service and are not redistributed by
Datagrok.

---

## 4. Docker container images (`dockerfiles/`)

Grokky ships two Docker images that run alongside the Datagrok platform.
These are built and distributed separately from the JavaScript bundle.
All listed components are under permissive licenses.

### `dockerfiles/claude-runtime` — Claude Agent runtime (Hono + WebSocket)

Wraps the official Claude Agent SDK and streams structured events to the
browser. Skills/knowledge-sync system included.

| Component                            | License    | Upstream                                          |
|--------------------------------------|------------|---------------------------------------------------|
| Node.js 20 (`node:20-slim` base)     | MIT (Node) | https://nodejs.org/                               |
| @anthropic-ai/claude-agent-sdk       | MIT        | https://github.com/anthropics/claude-agent-sdk   |
| @hono/node-server, @hono/node-ws     | MIT        | https://hono.dev/                                 |
| hono                                 | MIT        | https://github.com/honojs/hono                    |
| adm-zip                              | MIT        | https://github.com/cthackers/adm-zip              |
| yaml                                 | ISC        | https://github.com/eemeli/yaml                    |
| @zilliz/claude-context-mcp           | MIT        | https://github.com/zilliztech/claude-context      |

The image also clones `https://github.com/datagrok-ai/public.git` at build
time into `/workspace`; the cloned source is the same MIT-licensed `public/`
repository that this `CREDITS.md` lives in.

### `dockerfiles/mcp-server` — Datagrok MCP server (HTTP transport)

Exposes Datagrok operations (functions, files, projects, spaces, users) as
Model Context Protocol tools.

| Component                            | License    | Upstream                                          |
|--------------------------------------|------------|---------------------------------------------------|
| Node.js 20 (`node:20-slim` base)     | MIT (Node) | https://nodejs.org/                               |
| @modelcontextprotocol/sdk            | MIT        | https://github.com/modelcontextprotocol/typescript-sdk |
| @hono/node-server                    | MIT        | https://hono.dev/                                 |
| hono                                 | MIT        | https://github.com/honojs/hono                    |

---

## 5. Development-only dependencies

Tools used during the build/test cycle (not in the runtime tree, not bundled):
the `datagrok-tools` CLI transitively pulls in `puppeteer-screen-recorder`,
which references `@ffmpeg-installer/ffmpeg` (LGPL-2.1) and
`@ffmpeg-installer/win32-x64` (GPLv3). These binaries are **not** redistributed
as part of the published Grokky plugin and impose no obligation on users of
the plugin.

`@types/dom-speech-recognition` and `@types/dom-chromium-ai` are TypeScript
types-only (MIT) and not in the runtime bundle.
