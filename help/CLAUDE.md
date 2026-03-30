# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

This documentation serves multiple purposes:
1. Online help: https://datagrok.ai/help/datagrok/
2. In-application [Context Help](datagrok/navigation/panels/panels.md#context-help)
3. A resource for answering user questions using AI.

## Overview

This is the Datagrok documentation site, built with **Docusaurus 3**. Source files are GitHub-flavored Markdown (`.md`)
and MDX (`.mdx`). The Docusaurus config lives in `../docusaurus/`, not here.

Broken links are fatal: Docusaurus is configured with `onBrokenLinks: 'throw'` and `onBrokenMarkdownLinks: 'throw'`.

## Building and previewing

```bash
cd ../docusaurus
npm install
npm start             # Local dev server with hot reload
npm run build         # Full production build (catches broken links)
```

## Directory layout

| Directory      | Purpose                                         |
|----------------|-------------------------------------------------|
| `access/`      | Data access (databases, files, connectors)      |
| `collaborate/` | Collaboration features                          |
| `compute/`     | Scripting, compute environments                 |
| `datagrok/`    | Core platform concepts, navigation, solutions   |
| `deploy/`      | Deployment guides (Docker, AWS, K8s)            |
| `develop/`     | Developer docs, help-page style guides, how-tos |
| `explore/`     | Data exploration, multivariate analysis         |
| `govern/`      | Governance, access control, audit               |
| `learn/`       | Learning resources                              |
| `transform/`   | Data transformations                            |
| `visualize/`   | Viewers and visualization                       |
| `_access/`     | Internal tutorial content (not published)       |
| `_internal/`   | Selenium tests, admin docs, Dart utilities      |

## Sidebar ordering

- **Documents**: `sidebar_position` and `sidebar_label` in YAML front matter.
- **Folders**: `_category_.yml` file with `position` and `label` fields.

## Markdown conventions

### Front matter (required)

```yaml
---
title: "Short title"
sidebar_position: 2
---
```

Optional fields: `sidebar_label`, `description`, `keywords` (array), `format: mdx` (only when using
JSX), `unlisted: true`.

### Headers

- Never use H1 (`#`) in the body — `title` front matter renders as H1.
- H2 and H3 appear in the table of contents; keep them to 1-3 words.
- Sentence case only (capitalize first word and proper nouns).

### Links

All links must be **relative** to the markdown file. Never use absolute wiki URLs.

```markdown
[section in this doc](#section-header)
[doc in same folder](other-doc.md)
[doc in subfolder](subfolder/other-doc.md)
[doc in parent](../parent/other-doc.md)
```

### Images

Place in an `img/` subdirectory relative to the markdown file. Capture at 800×500px browser viewport. PNG and GIF
preferred.

```markdown
![alt text](img/screenshot.png)
```

### Content blocks

- Admonitions: `:::note`, `:::tip`, `:::warning`, `:::caution`, `:::important` — use sparingly.
- Collapsible sections: `<details><summary>Title</summary>...</details>`.
- Tabs: Import `Tabs` and `TabItem` from `@theme/Tabs` — only use when content varies by context (OS, language).
- Code blocks: specify language for syntax highlighting; optional `title="..."`.

### File naming

Lowercase, no spaces, kebab-case, short.
When a `.md` filename matches its parent folder name, clicking the folder shows that page.

## Documentation principles

- **Every page is page one.** Each page is self-contained — a reader can land anywhere and understand
  the topic without having read anything else. Cross-link rather than repeat.
- **Single source of truth.** The wiki is the authoritative reference for deploying, using, and
  troubleshooting Datagrok. If information doesn't exist here, add it — don't answer in Slack alone.
- **Live content.** Pages are automatically updated from platform core and plugin sources where
  possible. Tutorials are coded into the platform, then linked from docs — not duplicated as prose.
- **Link, don't repeat.** Reference cross-cutting concepts (functions, entities, viewers) by linking
  to their canonical page. Re-explaining creates drift.
- **User-oriented structure.** Organized by what users need: capabilities that build on each other →
  domain support (chem, bio, NLP) → fit-for-purpose apps.
- **Dual rendering.** The same markdown renders as the online wiki (datagrok.ai/help) and as
  in-app Context Help. Avoid complex HTML/MDX that won't work in the in-app widget.
- **Reference pages use tables.** Lookup-oriented pages present information in table format for
  quick scanning.
- **How-to pages are for complex procedures only.** Simple tasks get a sentence or a `<details>`
  block — not a dedicated page.

### Writing style

- US English, active voice, second person ("you").
- Concise: avoid nominalizations, filler words, passive constructions.
- Goal → Location → Action clause order for instructions.
- Bold UI elements: **Top Menu**, **Sidebar**, **Context Panel**.
- Oxford commas. No semicolons (use two sentences).
- Spell out zero through nine; use numerals for 10+.
- Punctuation outside quotation marks (British style); everything else US English.
- Contractions OK except noun+verb forms.

## Pre-commit checklist

1. Spelling and grammar (Grammarly or similar).
2. Front matter present with `title`.
3. No H1 in body.
4. No broken links or images.
5. Code snippets render correctly.
