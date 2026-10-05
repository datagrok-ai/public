---
paths:
  - help/**
  - docusaurus/**
  - docusaurus-static/**
---

# Documentation

`help/` is a submodule of datagrok-ai/help: commit and push doc changes in that repo, never the
`help` pointer in public (a bot bumps it on every help merge).

Documentation uses Docusaurus. Files are Markdown with YAML frontmatter:

```markdown
---
title: Page Title
sidebar_label: Short Label
sidebar_position: 5
---
```

Sidebar structure is defined in `docusaurus/sidebars.js`.
Images go in `docusaurus-static/images/`.
Internal links use relative paths without `.md` extension.
