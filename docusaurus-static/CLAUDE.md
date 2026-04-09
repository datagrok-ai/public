# docusaurus-static

Docusaurus site that serves Datagrok's documentation and landing pages. Embedded in the platform — served at `/` and `/help/` inside a running Datagrok instance.

## Purpose

- **`/help/*`** — documentation pulled from `../help/` (the `public/help/` markdown tree)
- **`/`** — custom landing/marketing pages (Hero, Body, Careers, Team, etc.)
- **Plugin release pages** — `generatePluginsPages.js` generates `plugins.json` from npm or local packages

## Key Directories

| Path               | Contents                                          |
|--------------------|---------------------------------------------------|
| `src/pages/`       | React pages (`index.tsx` = home, `new/`, `company/`) |
| `src/components/`  | Shared React components (home, card, team, careers) |
| `src/theme/`       | Docusaurus theme overrides (`MDXComponents.js`)   |
| `src/css/`         | Global CSS (`custom.css` — no navbar/sidebar in embedded mode) |
| `static/`          | Static assets (fonts, images, favicon)            |
| `../help/`         | Markdown docs source (routed as `/help/`)         |

## Scripts

```bash
npm run start    # generatePluginsPages.js create + docusaurus start
npm run build    # generatePluginsPages.js create + docusaurus build
npm run lint     # markdownlint-cli2-fix on help/**/*.{md,mdx}
```

## generatePluginsPages.js

Generates `../help/deploy/releases/plugins/plugins.json` and copies plugin CHANGELOGs into that folder.

- `node generatePluginsPages.js create` — creates `plugins.json` if missing (used at dev/build start)
- `node generatePluginsPages.js latest` — reads local `../packages/*/package.json` for versions
- `node generatePluginsPages.js` (default) — fetches all `@datagrok` packages from npm registry

## docusaurus.config.js Notes

- `onBrokenLinks/Anchors/MarkdownLinks` are all set to `'throw'` — broken links fail the build
- Docs root is `../help`, base path is `/help`
- Files matching `**/_*` are excluded from docs
- Color mode switch is disabled (light-only)
- Sidebar is hidden via CSS (embedded use case — sidebar navigation is external to this site)
