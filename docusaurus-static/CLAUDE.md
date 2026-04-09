# Docusaurus-Static (Simplified Documentation Build)

A stripped-down Docusaurus configuration, likely used for wiki or standalone doc builds separate from the main docusaurus site.

## datagrok.ai Site Architecture

The datagrok.ai website consists of two independently deployed parts on a single EC2 instance behind nginx:

1. **Landing** (Dart 1.x SPA) — marketing/corporate pages at `/`, `/company`, `/login`, `/solutions`, etc.
    - Source: `landing/` (private repo root)
    - Deployed by Jenkins → `/home/grok/landing/`
2. **Docusaurus** (`public/docusaurus/`) — documentation and API reference at `/help/*`, `/api/*`
    - Deployed by GitHub Actions → `/home/grok/docusaurus/`

Nginx routes decide which directory serves each URL path. The two deploys never overwrite each other.

For the full architecture, nginx routing rules, search (Typesense), CI/CD pipelines, and deployment details, see `landing/LANDING.md`.

## Key Files

| File                      | Purpose                                          |
|---------------------------|--------------------------------------------------|
| `docusaurus.config.js`    | Main config                                      |
| `sidebar-empty.js`        | Empty sidebar (no auto-generated sidebar)        |
| `generatePluginsPages.js` | Generates plugin release pages from CHANGELOG.md |

## Relation to Main Docusaurus

This is a lighter variant of `public/docusaurus/`. It uses a separate config and an empty sidebar, suggesting it serves a subset of the documentation or is used for a different build target (e.g., internal wiki). See the main `public/docusaurus/CLAUDE.md` for the full documentation site setup.
