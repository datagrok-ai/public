# Plates — BDD features

What Plates offers a user — the Plates app and its plate files — as features on
`@datagrok-libraries/bdd`. Plates has no help page; each feature says in its description what it
claims, from which source, and what it left out. Every step is the library's: this project has no
bindings of its own. It loads the `viewers` tier (grid and heat map readings, `no errors should have
been logged`).

| Feature | What it walks |
|---|---|
| `app/plates-app.feature` | Apps > Plates: the 40 demo plates with their properties as columns; Search Wells swapping its well inputs with the template; Search Plates narrowed by a number typed into a plate property; Templates with EDIT and CREATE NEW (not saved); Create Plate's analysis tabs and CREATE warning with no plate imported |
| `files/plate-text-file.feature` | a BMG PHERAstar .txt (`fixtures/bdd-bmg-plate.txt`) opened from the computer and from Browse > Files: the Summary heat map of a 96-well plate and each channel's values in its tab |

`features/*/MISSING.md` lists what these could not say yet: the writes into the plates database
(CREATE, template SAVE / Clone / Delete), well hover and click, Excel plates and the Assay Plates demo,
and the open questions.

Covered elsewhere: Search Plates and Search Analyses opening with their inputs —
`UsageAnalysis/bdd/features/browse/browse-apps-matrix.feature`.

## Running

```bash
cd public/packages/Plates/bdd
node ../../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=http://localhost:8888 node ../../../libraries/bdd/bin/grok-bdd.js run
```

The stand needs Plates published with its `plts` database; the package seeds its 40 demo plates on
first use. The file feature writes its fixture into My files and deletes it at feature end; the app
feature puts nothing on the server.
