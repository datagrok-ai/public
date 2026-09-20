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
| `panel/peptides-pane` | Renderer/metadata, real activity preview values, WebLogo selection |
| `sar/from-panel` | Pane launch, clustering settings, invariant map and distribution grouping |
| `sar/from-top-menu` | Dialog defaults, clustering settings, repeated optional-viewer lifecycle |
| `sar/weblogo-selection` | Exact source-row masks for click, Shift and Control/Command; dependent panes |
| `sar/tooltips` | Cell statistics, highlight cleanup, cluster statistics and selection |
| `sar/mutation-cliffs` | Independent pair counts, cell repaint, position chart, full export contents |
| `sar/similarity-threshold` | Fresh analyses at thresholds 10, 50, 75 and 90 |
| `sar/export` | Every invariant-map cell, every mutation pair/activity/delta, both source IDs |
| `sar/manual-alignment` | Apply, adjacent/end positions, Reset, selection against edited data |
| `sar/project-round-trip` | Non-default scaling/data/layout persistence and restored interactions |
| `entry/demo-dashboard` | Registered dashboard function, transformed activities and rendered viewers |
| `entry/landing` | Three demo buttons, their exact datasets/notations and side-panel state |

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
comparison. `HANDOFF.md` records the original survey and decisions; `docs/` contains the survey
and translation notes. Project persistence uses the public project API; it does not cover the
ribbon Save dialog. Dashboard invocation uses the registered function; gallery-card navigation
is separate. Dendrogram activation, per-cluster WebLogo glyph interaction, and Sequence Space
activation remain outside this translation's agreed scope.

## Run record

Every launch sets the MCL similarity threshold to 93 (the demo uses 94): MCL dominates the run
time and the default 70 spends minutes on this fixture to produce one cluster, while 93 clusters it
in seconds. The from-panel journey re-clusters at 90 and back at 93; the threshold outline runs
90, 93 and 96.

2026-09-15, macOS, stand `http://localhost:8889` (core and public at the working-tree state of
this translation): 15 tests (13 features; the threshold outline expands to 3) green serially
(53 s) and on four workers; before the threshold change the same suite took 3.6 min serially and
was green twice serially, twice on four workers and once headed. The `Widget: status providers`
and `Viewer: rendering` ApiTests categories pass on the same stand.
