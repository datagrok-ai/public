# Vendored code

Every directory under `vendor/` is copied from an upstream open-source project and frozen at
the recorded version. We own the code from the moment it lands here: adapt freely, never
chase upstream. Update this table with every addition or local change.

| Path | Upstream project | Version | License | Local changes |
|------|------------------|---------|---------|---------------|
| `signals-core/` | [@preact/signals-core](https://github.com/preactjs/signals) | 1.14.4 | MIT | none (verbatim `src/index.ts`) |
| `fontawesome-pro/` | [Font Awesome Pro](https://fontawesome.com) | 5.15.4 | Font Awesome Pro Commercial License (held by Datagrok) | `webfonts/fa-{light-300,regular-400,solid-900,brands-400}.woff2` copied verbatim from the platform's own assets (`core/client/xamgle/web/font/font-awesome/webfonts`); the same webfonts already ship publicly via `public/docusaurus/static/font/font-awesome`. The `.u2-standalone`-gated name→codepoint map in `../css/icons.css` is extracted from the platform's `fontawesome.css` |

## Ported (adapted, not vendored verbatim)

| File | Upstream | License | Notes |
|------|----------|---------|-------|
| `../src/components/list.ts` | [microsoft/vscode](https://github.com/microsoft/vscode) `src/vs/base/browser/ui/list/{listView,rangeMap}.ts` @ `5616258b8699` | MIT | render-range diff + row recycling ported onto u2 core; fixed height, no DnD/touch — details in the file header |
