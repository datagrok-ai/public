---
paths:
  - help/**
  - docusaurus/**
  - docusaurus-static/**
---

## Documentation

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
