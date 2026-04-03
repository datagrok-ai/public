---
paths:
  - libraries/**
---

## Library Development

Libraries are published under `@datagrok-libraries` scope. Each library has its own `package.json` and builds independently:

```bash
cd libraries/<lib-name>
npm install
npm run build              # Build the library
npm run lint               # ESLint check
npm run lint-fix           # ESLint auto-fix
npm run build-all          # Build this library and all its dependencies in order
npm run link-all           # Link local datagrok-api and other @datagrok-libraries/*
```

Libraries are consumed by packages. When modifying a library, rebuild it and re-link dependent packages:

```bash
grok link              # From within the package directory
grok link --unlink     # Revert to npm versioned dependencies
```
