---
paths:
  - libraries/**
---

## Library Development

Libraries are published under `@datagrok-libraries` scope. Each library has its own `package.json` and is a member of
the `public/` pnpm workspace (see `packages/BUILD.MD`). Never run `npm install` or `grok link` inside one.

```bash
grok setup                 # Once per checkout, at the repository root
cd libraries/<lib-name>
grok build                 # Build the library and what it depends on (tsc emit to dist/ with declarations)
pnpm run lint              # ESLint check
```

Libraries are consumed by packages through `workspace:^`. After modifying a library, run `grok build` in the
dependent package; Turborepo rebuilds the library first. Nothing to link.
