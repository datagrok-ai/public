# Docusaurus (Main Documentation Site)

Docusaurus 3.10 (React/Node.js) static site that serves the datagrok.ai documentation:
- `/help/*` — Platform documentation (Markdown from `public/help/`)
- `/api/js/*` — JavaScript API reference (TypeDoc-generated)
- `/api/py/*` — Python API reference (custom AST-based generator)

## datagrok.ai Site Architecture

The datagrok.ai website consists of two independently deployed parts on a single EC2 instance behind nginx:

1. **Landing** (Dart 1.x SPA) — marketing/corporate pages at `/`, `/company`, `/login`, `/solutions`, etc.
   - Source: `landing/` (private repo root)
   - Deployed by Jenkins → `/home/grok/landing/`
2. **Docusaurus** (this directory) — documentation and API reference at `/help/*`, `/api/*`
   - Source: `public/docusaurus/`
   - Deployed by GitHub Actions → `/home/grok/docusaurus/`

Nginx routes decide which directory serves each URL path. The two deploys never overwrite each other.

For the full architecture, nginx routing rules, search (Typesense), CI/CD pipelines, and deployment details, see `../landing/CLAUDE.md`.

## Key Files

| File                       | Purpose                                              |
|----------------------------|------------------------------------------------------|
| `docusaurus.config.js`     | Main config: plugins, theme, search, navbar          |
| `sidebars.js`              | Help sidebar (auto-generated from filesystem)        |
| `typedoc-sidebar.js`       | API sidebar (JS + Python combined)                   |
| `generatePluginsPages.js`  | Generates plugin release pages from CHANGELOG.md     |
| `generate_docs_py.sh`      | Runs Python API doc generator                        |
| `preprocess-readme.js`     | Rewrites JS API readme links to absolute URLs        |
| `docsearch.json`           | Typesense scraper configuration                      |
| `redirects.maps.json`      | URL redirect mappings                                |
| `broken_links.sh`          | Post-build broken link checker                       |

## Running Locally

```bash
npm install
cd ../js-api && npm install && cd ../docusaurus
npm run start-windows   # Windows
npm run start           # Linux/macOS
npm run start-debug     # Debug mode (warn instead of throw on broken links)
```

## Building

```bash
npm run build-windows   # Windows (sets NODE_OPTIONS for memory)
npm run build           # Linux/macOS
```

Build runs three pre-build scripts automatically: `generatePluginsPages.js`, `generate_docs_py.sh`, `preprocess-readme.js`.

## Deployment

Push to `master` in the `public` repo with changes in `help/`, `docusaurus/`, `js-api/`, or `python-api/`. GitHub Actions workflow (`.github/workflows/docusaurus.yaml`) handles build, deploy (rsync), and Typesense search re-indexing.
