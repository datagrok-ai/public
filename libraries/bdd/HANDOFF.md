# Working on the BDD tests — a handoff

For whoever picks this up next. `README.md` beside this file is the authoring reference (phrases,
steps, tiers, how a feature runs); `CLAUDE.md` is the engine's invariants; the `/bdd-translate` skill
is the translation procedure. **This file is the working guide**: the loop you actually run, the
judgement calls, the traps that each cost a full run to find, and where things stand.

---

## 1. Start here

A feature file is the test. It is written in the platform's own vocabulary, compiled into a
Playwright spec that the package commits, and the compiled spec is never edited by hand.

```bash
# one-time, per package, from the package directory
node ../../libraries/bdd/bin/grok-bdd.js init        # scaffolds bdd/
node ../../libraries/bdd/bin/grok-bdd.js link        # one @playwright/test per run

# the loop
node ../../libraries/bdd/bin/grok-bdd.js compile     # features -> generated/*.test.ts
node ../../libraries/bdd/bin/grok-bdd.js lint        # what every phrase resolved to
DATAGROK_URL=https://dev.datagrok.ai DATAGROK_SERVER=dev \
  node ../../libraries/bdd/bin/grok-bdd.js run --workers=3
```

`run` compiles with `--check` first, so a spec that drifted from its feature fails before a browser
starts. A global install or `npm link` gives you plain `grok-bdd` instead of the long path.

**The one principle, and it is not a slogan:** we test our own platform, not a black box. When a test
would wait, sleep, scan pixels or retry, a signal or a name is missing in the core — and the fix goes
there, not into the test. That is why the box plot suite runs in 14 s where the hand-written spec it
replaced took 49 s for one test. Everything in section 6 follows from this.

---

## 2. The map

| Where | What it is |
|---|---|
| `libraries/bdd/README.md` | authoring reference: element phrases, the step vocabulary, tiers, reading a failure |
| `libraries/bdd/CLAUDE.md` | the engine's invariants — read before changing the library itself |
| `libraries/bdd/bindings/common`, `bindings/platform` | the vocabulary every project gets: kinds, elements, gestures, the shell, the data steps |
| `libraries/bdd/bindings/tiers/viewers` | opt-in viewer vocabulary (properties, menus, hit areas, pixels, legend) |
| `libraries/bdd/src/runtime` | locate, gestures, assertions, the harness, the in-page viewer runtime |
| `<package>/bdd/features/*.feature` | the tests, in Gherkin |
| `<package>/bdd/bindings/*.ts` | that package's own vocabulary — only what is specific to it |
| `<package>/bdd/generated/*.test.ts` | compiled, committed, drift-gated; never hand-edited |
| `.claude/skills/bdd-translate/SKILL.md` | the procedure for turning an old spec into features |

Projects that exist today: `UsageAnalysis/bdd` (viewers + Spaces), `DiffStudio/bdd`, `Bio/bdd`,
`PowerGrid/bdd`, `U2Demo/bdd`, and the library's own smoke project.

---

## 3. The daily loop, and what "green" means

1. Write or change a `.feature`.
2. `compile` — it tells you what every phrase resolved to, and fails on a phrase that resolves to
   nothing. A resolution that surprises you is a bug in the phrase, not in the engine.
3. `run` the one feature you touched (`run generated/<name>.test.ts`), then the whole project.
4. Run the whole project **twice** before calling it done. The second run is where interference
   between features shows up — a feature that passes alone and fails in company is telling you its
   claims depend on state it did not state.

Useful flags: `--workers=N`, `--trace on`, `--video on`, `--repeat-each=2` (the flake test),
`--grep "<scenario>"`.

Green means: every feature of the project, twice, on the shared dev stand, with no feature left
depending on what an earlier one happened to leave behind.

---

## 4. Translating an old spec

The procedure is in `/bdd-translate`. What matters most in it:

- **Translate the case, not the old code.** The TypeScript spec's workarounds are not requirements:
  a canvas hash becomes the viewer's own "should have repainted", a pixel scan becomes a hit area, a
  sleep becomes a signal. If the old spec worked around something, find what it was working around.
- **Every scenario ends on `no errors should have been logged`.** A scenario owns its error floor.
- **Then review each feature backwards**: read the feature against the spec it replaced, line by
  line, and ask of every claim "what would have to break for this to fail?". A claim that cannot fail
  is worse than no claim. This round is where most of the real findings came from — a per-area repaint
  claim that measured the whole canvas, a tooltip reader that matched the previous hover's text, a
  negative tree claim that was true of a closed group whatever the server held.
- **State what you did not translate, and why**, in the feature's own description. A reader must not
  have to diff against the old spec to learn that.

---

## 5. Where a binding goes, and how to find a phrase

**The rule:** reuse the existing vocabulary first. A binding that could serve anywhere else belongs in
`libraries/bdd`, with the full Given/When/Then family and a scoped phrase
(`the {string} tree node inside of browse tree is expanded`, not a bare `tree node`). Only what is
genuinely specific to one package stays in that package's bindings — its app's routes, its own
widgets, its own server artefacts.

**Finding the phrase** — never guess a selector, probe the page. The recipe that pays for itself:

```js
// scratchpad/probe.mjs — createRequire so it uses the library's Playwright
const require = createRequire('<...>/libraries/bdd/package.json');
const {chromium} = require('@playwright/test');
const page = await (await browser.newContext({storageState: '<pkg>/bdd/e2e/.auth.json'})).newPage();
// goto, then dump what the engine will see: tag, classes, data-u2*, role, aria-*, name, value, label
```

Dump `[name^="input-host-"]`, `.d4-ribbon-item`, `[role="treeitem"]`, the open popup — whatever the
phrases will name. Write probes with the Write tool, not a heredoc (see the traps table).

Phrase resolution order: a registered whole phrase → an ordinal → a split at the first scope word
(`X in|inside|within|on|of Y`) → a registered element → a generic kind by suffix. Dart names come
from `dartNames` templates (`input-host-{q}`, `tree-{q}`, `button-{q}`, `div-section--{q}`), with the
qualifier's spaces as dashes.

---

## 6. When it fails

A failure tells you the feature line, the step, one sentence, and what was on the page instead
(`visible inputs: …`, `it has: …`). Read that before opening a trace; it is usually enough.

- **A throw inside `expect.poll` ends the poll.** A read that may hit a transient state returns
  `false` and keeps the reason for the message — it must not throw.
- **"Before" is before the last change**, never after the last check: the than-before family does not
  move the snapshot when it passes.
- **Instrument before blaming the product.** Every single time a wall looked like a platform bug, it
  was a wrong selector, a race, a collapsed group, or my own wrong assertion — five false "product
  bugs" before that sank in. Measure first: a probe that prints the state, a `--repeat-each=2` run, the
  same feature alone versus in company, the same suite on `origin/master` versus on your branch.
- **A/B against master is the cheapest attribution you have.** `git checkout origin/master -- libraries/bdd`,
  rebuild, run, compare, restore. That is how the Spaces failures were shown to be master's own and
  the nine viewer failures to be the stand's.
- **Never file a ticket from a run.** Propose the finding with its evidence; Olesia files it after
  walking the steps by hand.

---

## 7. Traps that each cost a run (all measured on dev)

| Trap | What to do |
|---|---|
| A compute form's parameter switch is **not inside** the input: a child of the host in the sensitivity view, the **preceding sibling** in the fitting view, and its own host has no name | use `user switches on {element}` — it searches itself, inside, then back over the siblings |
| Switching a parameter on **replaces** its input: `FFox` becomes `FFox min`/`FFox max`, and the switch moves to the min | claim the min/max, not the input you just switched |
| A sensitivity run with no parameter on ends with **one** viewer and says nothing | switch parameters on first |
| The fitting Run icon stays grey until every varied parameter has **both bounds**; the reason is only in its tooltip. FKox, which the manual case names, carries none in the Bioreactor model | vary FFox (0.15–0.25), or fill the bounds |
| The Model Hub gallery is on the page, visible and **empty for ~6 s** after `Compute2:modelCatalog` returns | wait for cards, never for the element |
| The platform's script view is **CodeMirror 5** (`.CodeMirror`), and running the script replaces the editor with a RichFunctionView | tag and save **before** Run |
| The browse tree's Spaces group is closed unless something opened it, and the product's reveal of a new space is unreliable under load. A closed group made every `should be absent` claim **vacuous** | use the claims that open the group themselves |
| A dialog that commits to the server closes when the server answers — 6–18 s, straddling the shared 15 s budget | `the {string} dialog should close` |
| A Dart `SwitchInput` hides its real checkbox and publishes state only as `ui-input-switch-on`; a ribbon button says "disabled" with a class | read what exists, and consider adding the aria state in core |
| A viewer that answers `getWidgetStatus` with **no readings** is an old client, not a broken test | check the stand's build before chasing it |
| A package's webpack type-checks its `bdd/` folder unless tsconfig excludes it | `"exclude": [..., "bdd"]` |
| The library needs **Node 20** (its bindings import Playwright at load); the libraries workflow pins 18 unless `engines.node` says otherwise | keep `engines` in the library's package.json |
| The Bash tool collapses `\\` inside heredocs, and `cd` persists between calls | write code with Write/Edit; use absolute paths or `git -C` |

---

## 8. Discipline

- **Branch only.** BDD work is committed and pushed to its own branch so colleagues can see it;
  merging to master needs an explicit decision. Stage by explicit path — the tree carries other
  sessions' work.
- **Clean up what a test creates.** Every step that writes to the stand registers its undo with
  `atFeatureEnd` and deletes exactly what it added. Two earlier suites left 23 `.ivp` files and three
  model-tagged scripts on dev, which is how `PK-PD(22).ivp` came to exist.
- **A sizable change gets a one-line `## v.next` entry** in the package's `CHANGELOG.md`.
- **`@known-failure`** marks a scenario that states a defect the product has: its failure does not
  fail the test and its **passing** does, so the tag must go when the bug is fixed. Never use
  `test.fail()` and never soften a claim to stay green.
- Before the pull request: the library's unit tests (`npm test` in `libraries/bdd`), `lint`,
  `compile --check` in **every** project (a library change can break another project's phrases), and
  the suites of the areas you touched.

---

## 9. Where things stand (11 Sep 2026)

| Area | State |
|---|---|
| Viewers (`UsageAnalysis/bdd/features/viewers`) | 121 features. 9 fail on dev for the stand's own reason: the client it serves (1.28.0) reports no readings for the word cloud and the map, and the Forms viewer comes from the PowerGrid build on the stand. They need dev rebuilt with current core and PowerGrid republished |
| Spaces (`UsageAnalysis/bdd/features/spaces`) | 10 features, green. Merged to master, then stabilised again (the tree-group claims above) |
| Diff Studio (`DiffStudio/bdd`) | 9 features, 40 scenarios, green. Covers all eight TestTrack specs, including the three cases the old suite left to a human |
| Bio, PowerGrid, U2Demo | projects exist from earlier rounds |
| Open pull request | [public#4073](https://github.com/datagrok-ai/public/pull/4073) — Diff Studio plus the Spaces claims and five CI fixes |

Not translated yet, in rough order of value: Browse/Demo apps, PowerPack, Chem, EDA, Peptides,
Tutorials, Users/Groups/Roles, StickyMeta, Home widgets, the Apps matrix.

Two known-deliberate gaps in Diff Studio, both stated in the features: a second browser tab (a
feature gets one page, so the address is reloaded instead) and the manual case's `FKox 0→3` (the model
gives FKox no bounds, so the fit cannot run on it).

---

## 10. Before you call an area done

- Every scenario of the area's features passes, twice, on the dev stand.
- Each feature passes **alone** and **in company** (run the whole project with 4 workers).
- Every claim would fail if the thing it names broke — you checked the negative ones especially.
- Nothing the suite created is left on the stand.
- The features say what they did not translate, and why.
- `compile --check` is clean in every project, the library's unit tests pass, `lint` is clean.
- The area's md specs, if they carried a "manual only" list, have been re-read: most of those items
  are automatable once you stop reproducing the old spec's workaround.
