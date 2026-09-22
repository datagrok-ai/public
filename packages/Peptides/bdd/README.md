# Peptides BDD tests

The suite translates the package's hand-written Playwright tests into Gherkin. Features use the
shared `@datagrok-libraries/bdd` vocabulary and real mouse/keyboard interactions with named UI
elements. Package bindings add peptide-specific arithmetic checks and entry-point setup.

## Run

From `public/packages/Peptides`:

```sh
npm install
npx grok-bdd link
npx grok-bdd compile
npx grok-bdd run --workers=1
npx grok-bdd run generated/sar/tooltips --workers=1 --headed
```

`grok-bdd link` makes the package use the BDD library's Playwright installation. Run it after
installing dependencies if Playwright reports that `@playwright/test` was required a second time.
Features and bindings are the sources; regenerate `bdd/generated/` instead of editing specs.

The runner defaults to `http://localhost:8888`; point it at another stand with `DATAGROK_URL`
(`DATAGROK_URL=http://localhost:8889 npx grok-bdd run --workers=1`). It signs in with the dev key
of the matching `~/.grok/config.yaml` server, else through the login form with `DATAGROK_LOGIN` /
`DATAGROK_PASSWORD` (`admin`/`admin` by default). Authentication state is local, not checked in.

The stand needs Peptides and its Bio/EDA dependencies published. A change to MCL in
`libraries/ml` reaches the stand only through **EDA**, which registers that viewer; publishing
Peptides alone does not replace EDA's copy. From each package directory: `grok build` and
`grok publish localhost`.

The served core and JS API must include `DG.Widget.addStatusProvider`, which Peptides uses to
publish the WebLogo glyphs it draws in the native grid's headers as hit areas.

## Coverage

| Feature | Main checks |
|---|---|
| `panel/peptides-pane` | Renderer/metadata, real activity preview values, WebLogo selection, Clusters vs Generate clusters, Manual Alignment without an analysis |
| `sar/from-panel` | Pane launch, completed MCL threshold, clustering settings, invariant map and distribution grouping, cliff glyphs |
| `sar/from-top-menu` | Dialog defaults, dock layout, clustering settings and their effect on the clusters, settings-dialog state, repeated optional-viewer lifecycle, Sequence space, Dendrogram on and off |
| `sar/weblogo-selection` | Exact source-row masks for click, Shift and Control/Command; Distribution and Selection panes after each |
| `sar/tooltips` | Cell statistics, cliff glyphs, highlight cleanup, cluster statistics and selection |
| `sar/mutation-cliffs` | Independent pair counts, p-values, a cliff cell drawn and a cliff-free cell left blank, position chart, full export contents |
| `sar/similarity-threshold` | Fresh analyses at thresholds 10, 50, 75, 90, 93 and 96 on 200 peptides and at 90 on all 647; header and map selection after each |
| `sar/export` | Every invariant-map cell, every mutation pair/activity/delta, both source IDs |
| `sar/manual-alignment` | Apply, adjacent/end positions, recomputed statistics, Reset, selection and panes against edited data |
| `sar/project-round-trip` | Non-default scaling/data/layout/selection persistence and restored interactions and panes |
| `sar/default-launch` | The dialog defaults (threshold 70, inflation 1.4, no scaling) run once end to end |
| `entry/demo-dashboard` | Registered dashboard function, transformed activities, WebLogo headers and rendered viewers |
| `entry/landing` | Exactly three demo buttons, their exact datasets/notations and side-panel state |

The source fixture is `System:DemoFiles/bio/peptides.csv`: 647 rows, 17 separator positions,
22 non-gap monomers and 6253 unique single-mutation row pairs. Export oracles derive expectations
from the original sequence strings and activities, independently of the viewer's caches.

Journeys share their Background and run scenarios in order. Each scenario checks its own console
and balloon errors because the harness clears those between scenarios. `peptides-sar-ready`
signals completed launch/settings work; viewer rendering is synchronized through actual pending
work and render events. DOM clicks on radio buttons need an explicit viewer snapshot before a
repaint comparison.

## Boundaries and review

The original `playwright/` specs and `public/playwright-public/Peptides/` descriptions remain for
comparison; `docs/` contains the original survey, decisions and translation notes. Each feature's
description states what it leaves out and why. Project persistence uses the public project API;
it does not cover the ribbon Save dialog. Dashboard invocation uses the registered function;
gallery-card navigation is the Browse suite's. Per-cluster WebLogo glyph interaction is not
exposed as hit areas (the Logo Summary Table reports whole cells only).

## Run record

Every launch but one sets the MCL similarity threshold to 93 (the demo uses 94): MCL dominates the
run time and the default 70 spends a long time on this fixture, while 93 clusters it in seconds.
`sar/default-launch` runs the defaults once. The from-panel journey re-clusters at 90 and back at
93; the threshold outline runs 10, 50, 75, 90, 93 and 96 on 200 peptides, plus 90 on all 647.

Known failures, each a candidate finding waiting for its ticket (stated above the scenario):
Sequence space checked in the settings embeds nothing; the settings reopen with Dendrogram
unchecked while the tree is shown; unchecking Dendrogram leaves the tree.

2026-09-22, dev.datagrok.ai (core master, Peptides and EDA published from this checkout as debug
builds): 19 tests (13 features; the threshold outline expands to 6) green three times in a row on
three workers, 1.9–2.1 min each. The JSON reports show each known failure failing at its own step.
