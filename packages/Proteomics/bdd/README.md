# Proteomics — BDD features

The Proteomics package's import as features on `@datagrok-libraries/bdd`. Each feature says in its description
what it claims and what it left out. Every step is the library's: this project has no bindings of its own; it
loads the `viewers` tier for the error floor (`no errors should have been logged`).

| Feature | What it walks |
|---|---|
| `import/import-files.feature` | the left-bar main menu's Proteomics > Import: MaxQuant... on `files/demo/proteinGroups.txt` (contaminant and reverse hits left out, protein IDs typed), Generic Matrix... and its column-mapping dialog |

The file is chosen in the browser's file chooser the menu opens, straight from the package's `files/demo/`
(the step resolves the path against `bdd/`). `features/import/MISSING.md` lists the questions and what is not
written yet.

## Run

From `public/packages/Proteomics/bdd`:

```sh
node ../../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=http://localhost:8888 node ../../../libraries/bdd/bin/grok-bdd.js run --workers=2
```
