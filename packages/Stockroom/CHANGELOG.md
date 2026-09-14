# Stockroom changelog

## v.next

* GROK-20753: Introduced Stockroom — a chemical stockroom on the GHS classification and the zero-code reference app for entity-mapped domain schemas: twelve tables declared in `databases/stockroom/schema.json` (constraints, ref filters, searchable columns, custom permissions, a self-referencing location tree, N:N hazards, a file column, auto-numbered labels), the UNECE GHS vocabulary (29 hazard classes, 80 H-statements, 97 P-statements) plus a small demo stockroom as seed scripts, a three-line `package.ts` over `domains.table(...).app()`, a `dg-ui/1` spec for the configuration tier, and a schema smoke test
