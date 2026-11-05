# Compute2 — BDD features

The Model Hub's own UI as features on `@datagrok-libraries/bdd`. Each feature says in its description what it
claims and what it left out. Every step is the library's: this project has no bindings of its own; it loads the
`viewers` tier for the error floor (`no errors should have been logged`).

| Feature | What it walks |
|---|---|
| `model-hub/model-hub-tools.feature` | the Model Hub opened from Apps > Compute, its gallery search narrowing the catalog, Tools > Compare Runs... opening an empty Run Comparison view |

`features/model-hub/MISSING.md` lists Browse Runs... (not on the local stand yet) and the run comparison
scenario that needs saved runs and a named model picker.

Covered elsewhere: opening the Model Hub from Browse and a saved model's preview and menu
(`UsageAnalysis/bdd/features/browse/browse-apps-and-dashboards.feature`), a Diff Studio model run from the
catalog (`DiffStudio/bdd/features/catalog.feature`), the compute tutorials (`Tutorials/bdd/features/compute/`).

## Run

From `public/packages/Compute2/bdd`:

```sh
node ../../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=http://localhost:8888 node ../../../libraries/bdd/bin/grok-bdd.js run --workers=2
```
