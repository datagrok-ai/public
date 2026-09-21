---
feature: connections
realizes_atlas: []   # untyped: environment-hygiene check, not a product-feature type
realizes: [views.databases]
priority: p2
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  Environment-inventory judgment on the public stand: which demo databases and
  packages are deployed is a property of that environment, not of the product,
  so a stand-agnostic automated spec cannot assert it.
related_bugs: []
---

# Databases node — public environment hygiene

Verifies that the public environment shows only the intended demo databases
and that no test packages or test connections leak into the user experience.
This scenario is public-environment-specific by design.

## Preconditions

- Run on the public environment; Databases is visible in Browse.
- Expected demo databases: **Chembl, Northwind, Starbucks, World**, and
  optionally **SureChEMBL**.

## Steps

1. Go to **Browse > Platform > Packages** and search for `DbTests`, then
   `ApiTests` — neither package is installed
2. Go to **Browse > Databases** and expand the node for the first time —
   only providers with at least one active connection are listed; providers
   without any connection do not appear
3. Verify the visible demo databases are exactly: Chembl, Northwind,
   Starbucks, World — plus SureChEMBL if deployed, named exactly
   "SureChEMBL"
4. Verify **Athena** is absent, and no listed connection shows an error or
   disconnected state when expanded

## Cleanup

Read-only scenario — nothing to clean up.

---
{
  "order": 10,
  "datasets": []
}
