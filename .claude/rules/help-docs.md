---
paths:
  - docusaurus/**
  - docusaurus-static/**
---

# Documentation

Help pages live in the datagrok-ai/help repository, a sibling of `public/` in the datagrok umbrella
(`../help`): commit and push doc changes there. A help merge to master redeploys the site
(help's `bump-public` workflow dispatches `docusaurus.yaml` here).

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
