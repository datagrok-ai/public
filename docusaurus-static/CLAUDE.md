# Docusaurus-Static — OBSOLETE

> **This directory is obsolete.** All work should be done in `public/docusaurus/` instead.
> `docusaurus-static/` was previously used for experimental landing page work (Team page, company pages)
> but is not wired into any CI/CD pipeline and changes here are never deployed to datagrok.ai.

## Where to work instead

- **Team page, company pages, React components**: `public/docusaurus/src/`
- **Documentation**: `public/docusaurus/` — deployed via GitHub Actions (`.github/workflows/docusaurus.yaml`) on push to `master`

See `public/docusaurus/CLAUDE.md` and `landing/LANDING.md` for the full architecture.
