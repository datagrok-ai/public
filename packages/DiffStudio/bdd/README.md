# Diff Studio — BDD features

The eight TestTrack cases of Diff Studio (`packages/UsageAnalysis/files/TestTrack/DiffStudio/*.md`)
and the Playwright specs beside them, as features on `@datagrok-libraries/bdd`. Each feature says in
its description what it translates and what it left out.

| Feature | Case | What it walks |
|---|---|---|
| `hub.feature` | open-model step 1 | the app's hub, a library card, the Open model icon's Library menu |
| `open-model.feature` | open-model | Bioreactor: tabs, the Facet colours, the slider, the address, Process mode's cascade |
| `scripting.feature` | scripting, M-1.5, M-1.6 | Edit, `</>`, the script's tabs and grid rows, save with `//tags: model`, the Model Hub, Refresh |
| `sensitivity-analysis.feature` | sensitivity-analysis | the Edit toggle both ways, the switches, a run over three parameters |
| `fitting.feature` | fitting | the switches and bounds, a target table from the stand, a fit whose loss never increases |
| `catalog.feature` | catalog | save to the library (claimed by the app's own event), the Model Hub, the model run from there |
| `cyclic-models.feature` | cyclic-models | PK-PD: the clickers, the tooltips |
| `stages.feature` | stages | Acid production: a stage duration, the tooltips |
| `files-and-sharing.feature` | files-and-sharing | pk.ivp from the Files tree, the slider, the address loaded again |

## Running

```bash
cd public/libraries/bdd && npm ci && npm run build && npx playwright install chromium
cd ../../packages/DiffStudio && npm ci && npx grok-bdd link
npx grok-bdd run --reporter=list            # the library's 4 workers; --workers=2 when the stand falls behind
npx grok-bdd run generated/fitting.test.ts  # one feature
```

Nine features, 42 scenarios: 1.9 min on three workers against a local stand (2026-09-11), twice
in a row. The stand needs `DiffStudio` and `Compute2` published, the library's files under
`System:AppData/DiffStudio/library` (the package's `files/library`), and a dev key for the login
(`libraries/bdd/README.md`, "What the stand needs"). Every feature leaves the stand as it found it:
the library file a save adds and the script the script view saves are deleted when the feature
ends.

## What the platform gave these features

- The model's charts are real viewers, in the model view and in the function views the Model Hub
  and the script view open — the library reaches a viewer outside a table view through
  `DG.Widget.find`, and a function view's tabs (`.dockspan-tab-handle`, in a shadow root) are
  `tab`s. A claim names the tab it reads: the viewers of the other tabs stay in the DOM with no
  rectangle.
- A compute form's parameter switch is a Dart `SwitchInput` with `role="switch"` and
  `aria-checked`, which is what `user switches on {element}` and `should be switched on` read.
- The Facet plot is small multiples, several line charts: `the canvases of open tableview should
  be painted in at least N colors` counts across them.
