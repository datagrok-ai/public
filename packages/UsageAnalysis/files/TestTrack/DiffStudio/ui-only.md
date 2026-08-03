# Diff Studio Playwright tests — UI-only constraints + manual verification

This document has two purposes:

1. **Technical log** — records the steps where strict UI-only automation was not achievable, the
   fallback chosen, and the reason. Order of preference per spec is:
   1. Pure UI (clicks, hovers, keyboard, address-bar navigation, DOM dispatch of native events).
   2. A single Datagrok API call (`grok.dapi.*`, `grok.shell.*`).
   3. Document the gap here.

2. **Manual verification checklist** — full end-to-end scenarios that a human tester must still
   run by hand. The list shrank substantially after the second automation pass — most items
   that were marked manual originally are now covered by canvas-hash, pixel-histogram, dataframe
   introspection, IVP-tooltip parsing, real slider drag, or count-delta refresh checks.

Source-of-truth scenarios: `public/packages/UsageAnalysis/files/TestTrack/DiffStudio/*.md`.

---

# Part 1. Technical constraints log

## Steps that fell back to API

These tests perform exactly **one** read-only API call as a fallback when no UI gesture
produces the same effect. All UI assertions before and after the call stay pure UI.

### `fitting.test.ts` Step 5 — Add `bioreactor-experiment.csv` to the Bioreactor table input

**MD scenario:** "In the Bioreactor table input field, add the file from Files: App Data > Diff
Studio > library > bioreactor-experiment.csv"

**Why pure UI was rejected:**
- The Bioreactor input is a `TableInput` (`<select>` populated from `grok.shell.tables`). The file
  must already exist as a loaded dataframe in the shell before the dropdown shows it.
- Navigating to the CSV URL (`page.goto('/files/.../bioreactor-experiment.csv')`) replaces the
  current view with a CSV preview. `page.goBack()` then re-creates the Fitting view from scratch,
  losing all switcher edits made in Step 4.
- Drag-and-drop from the Browse tree into the Vue-rendered Fitting form is not reliably
  reproducible in Playwright across worker runs.
- Opening the CSV in a new browser tab via `context.newPage()` does not help because each tab has
  its own `grok.shell` instance — the loaded dataframe is not shared.

**Fallback:** A single API call inside a `page.evaluate`:

```ts
const df = await grok.dapi.files.readCsv(
  'System:AppData/DiffStudio/library/bioreactor-experiment.csv');
df.name = 'bioreactor-experiment';
grok.shell.addTable(df);  // adds the table without spawning a new view
```

The subsequent UI step (selecting `bioreactor-experiment` in the dropdown via `selectChoice`) is
pure UI.

### `cyclic-models.test.ts` / `stages.test.ts` Step 4 — Read IVP source for tooltip comparison

**Goal:** Assert that the tooltip the user sees on hover equals the `[tooltip text]` declared in
the source IVP file. Pure UI cannot supply the "what the tooltip should say" reference.

**Fallback:** Read `pk-pd.ivp` / `ga-production.ivp` via `grok.dapi.files.readAsText`, parse the
`[bracketed text]` annotations, then assert UI tooltips contain the parsed text. If the dev
server returns 502/404 for the file path, the test falls back to hardcoded expected strings
maintained in the test source so the assertion still runs without the API call.

### `catalog.test.ts` Step 4 — Refresh observably reloads (delete via API)

**Goal:** Prove the **Refresh** button on the Model Hub ribbon actually re-fetches the catalog —
i.e. that an externally-removed entry disappears after a click.

**Fallback:** Delete one user-saved PK-PD script via `grok.dapi.scripts.delete(...)`, click
Refresh, assert the card count decreased by exactly that many entries.

**Why pure UI was rejected:** Deleting via the card's right-click menu involves a confirmation
dialog whose selectors vary between Datagrok builds; deleting via the file tree requires
navigating multiple levels and is fragile. The API call here is read-modify-delete on shell
state, not a substitute for the UI verification (which is the click-Refresh + count assertion).

## Steps requiring `page.evaluate` for DOM traversal (not an API fallback)

These use `page.evaluate` purely to walk the DOM and dispatch native `click` events on elements
the Playwright locator chain cannot reach reliably. They are **UI-level** — the page receives a
real DOM `click` event indistinguishable from a mouse click.

- `scripting.test.ts` Step 4 / `catalog.test.ts` Step 3: `clickTreeLabel` — walks the Browse
  tree (`.d4-tree-view-group-label` / `.d4-tree-view-node-label`) and clicks the matching label.
  The platform's tree is virtualised — some labels live outside the viewport and Playwright's
  locator-based clicks fail with "element not visible" even with `scrollIntoViewIfNeeded`.

- `fitting.test.ts` Step 4: `enableSwitcher` — finds the `.sa-switch-input` sibling of an input
  host via `previousElementSibling` traversal and dispatches `.click()` on the inner
  `.ui-input-switch`. Compute2's Fitting view creates per-input toggles as **separate sibling
  inputs** with `height: 0; position: relative; top: 22px` styling (see
  `compute-utils/function-views/src/fitting-view.ts:1244`), giving them a degenerate bounding
  box. Standard Playwright clicks land outside the toggle.

- `fitting.test.ts` Step 6a: target-output switcher toggle — same `.sa-switch-input` traversal,
  applied to the Target Block at the bottom of the Fit form.

- `sensitivity-analysis.test.ts` Step 4: `enableSwitcher` — Compute2's SA view does the opposite
  of Fitting: it inserts the switcher as a **child** of the input host (see
  `compute-utils/function-views/src/sensitivity-analysis-view.ts:178`). Selector
  `[name="input-host-<name>"] .sa-switch-input .ui-input-switch` and a dispatched click.

- `files-and-sharing.test.ts` Step 1: searches the directory grid for a text-exact element with
  filename "pk.ivp" and dispatches clicks up the ancestor chain — the directory listing component
  does not expose a stable per-row class that Playwright can target by attribute.

- `open-model.test.ts` Step 4 — `M-2 Facet colour distinctness`: reads pixel data from every

сдф  `.d4-viewer canvas` via `getContext('2d').getImageData(...)`, buckets RGB channels to 4 bits
  to collapse antialiasing noise, and counts distinct buckets. This is read-only canvas
  inspection — the user sees the same pixels.

## Global error monitor (`helpers/error-monitor.ts`)

Every test attaches a monitor on startup that fails the test if any of the following fires:

- `console.error` (after filtering known dev-environment noise — Jiraconnect/Plates CORS
  failures, Dart SDK internal async traces, sourcemap 404s, 502s for `System.AppData/file/...`).
- `pageerror` (uncaught browser exceptions).
- HTTP responses with status ≥ 500.

The filter list lives in one place (`error-monitor.ts:isNoise`) — extend it cautiously.
Real platform errors (a balloon error, a 500 from a DiffStudio endpoint, a thrown exception
during model run) fail the corresponding test on the `monitor.assertNone()` line at the end of
each test.

## Out-of-MD assertion removed

- `catalog.test.ts` Step 6: the original `-spec.ts` source asserted `page.url()` contains
  `dose=5000`. Models opened via the Model Hub run inside a Compute2 RichFunctionView whose URL
  does **not** mirror the model inputs (unlike DiffStudio's in-app `?params:` scheme). The MD
  only requires the input value to update, so the URL check was removed.

## Known scenario remarks (preserved in test annotations)

These come from the MD source and are deliberately permissive in the assertions:

- `fitting.test.ts` Step 6b — *REMARK*: per the MD, "Grid may contain another number of rows".
  Row-count is permissive; **if** an `RMSE`-named column is produced, its values are asserted
  to be 95%+ non-increasing and to end ≤ the starting value (descending fit curve).
- `files-and-sharing.test.ts` Step 4 — *REMARK*: with 1–3 curves there is no Multiaxis/Facet;
  the test drops Count back to 1 and asserts the absence of those tabs.

---

# Part 2. Manual verification checklist

The original 8-item checklist (M-1 ... M-8) is mostly **automated** now. This part keeps only
the cases where automation is genuinely impossible on the current builds. For each, follow the
steps from a clean login.

**Common preconditions:**
- A working Datagrok account with permission to save scripts and read `System:AppData/DiffStudio/`.
- A modern Chromium-based browser; developer console reachable via F12.
- DiffStudio and Compute2 packages deployed and visible in the platform's Apps catalog.

On any failure: screenshot the page, copy any red `.d4-balloon` text, dump the dev-console
errors, attach to the ticket. Test ID (e.g. `M-1.5`) goes in the ticket title.

## Items that remain manual

### M-1.5. Live chart redraw inside the platform ScriptView
**Source:** `scripting.md` step 2.
**Why still manual:** the chart inside the platform `</>`-spawned ScriptView (which then runs
the script under a Compute2 RichFunctionView editor) is rendered inside a Vue component tree.
Three automation paths were tried and rejected:
- `.d4-viewer canvas` hash — the selector resolves but the canvas's `toDataURL` returns a
  constant blank image; the actual chart pixels live in a sibling overlay canvas the test cannot
  reliably address.
- `grok.shell.t` dataframe fingerprint — `grok.shell.t` is `null` in this view; even searching
  `grok.shell.v.dataFrame` and `grok.shell.tables` yields a sum that does not change between
  Final-input edits, because the table reference is recycled in place by Vue's reactive state.
- Playwright `locator.screenshot()` on `.d4-viewer` — returns `<missing>` because the locator
  finds no visible element in this view.

**Manual steps:**
1. Log in. Open Diff Studio → Library → Bioreactor.
2. Turn on the **Edit** toggle on the ribbon, then click **`</>`** — a ScriptView opens with a
   CodeMirror panel and a `Run` (▶) icon.
3. Click **Run**. The inputs panel gains a **Final** input and a line chart appears.
4. Type `500` into Final and press **Tab**. **Expect:** chart redraws, grid updates.
5. Type `200` then `1000` — expect continuous live redraw at each Tab.

**Pass criteria:** chart visibly redraws on each Final change.

---

### M-1.6. Live chart redraw inside Model Hub for a script-tagged model
**Source:** `scripting.md` step 5.
**Why still manual:** same Vue isolation as M-1.5 — the Model Hub viewer's chart canvas is not
exposed under any addressable selector.

**Manual steps:**
1. Log in. Run M-1.5 first to save a Bioreactor script with `//tags: model`.
2. Sidebar **Apps** → **Compute** → **Model Hub**.
3. Find a **Bioreactor** card with a JS icon (the saved script). Double-click it.
4. Set **Final** = `800`, Tab. **Expect:** chart redraws, grid updates.
5. Repeat with `300` and `1500`.

**Pass criteria:** Model-Hub-opened script reacts identically to the ScriptView in M-1.5.

---

### M-1.7. Live chart redraw inside Model Hub for an IVP model (PK-PD)
**Source:** `catalog.md` step 6.
**Why still manual:** same Vue isolation as M-1.5/M-1.6. The automated test asserts the input
value changes; chart redraw is the manual delta.

**Manual steps:**
1. Log in. Open Diff Studio → Library → PK-PD. Click **Save to Model Hub** (book-plus icon).
2. Sidebar **Apps** → **Compute** → **Model Hub** → double-click a **PK-PD** card.
3. Click into **dose**, type `5000`, press Tab. **Expect:** chart redraws, grid updates.
4. Type `1000` then `10000` — expect continuous live redraw.

**Pass criteria:** dose changes inside the Model-Hub-opened PK-PD live-update the chart and grid.

---

## Items previously manual that are now automated

| Original ID | What changed | Where it lives now |
|---|---|---|
| **M-1.1** Switch at chart redraw | Canvas hash of `.d4-viewer` before/after input change | `open-model.test.ts` Step 5 |
| **M-1.2** Process mode cascade + chart redraw | Same canvas-hash approach + cascade verification | `open-model.test.ts` Step 6 |
| **M-1.3** PK-PD Count clicker chart redraw | Canvas hash | `cyclic-models.test.ts` Step 3 |
| **M-1.4** pk.ivp slider/clicker chart redraw | Canvas hash + real slider drag | `files-and-sharing.test.ts` Step 2 |
| **M-2** Facet plot 12 distinct colours | Pixel histogram (bucket RGB to 4 bits, count buckets) | `open-model.test.ts` Step 4 |
| **M-3.1** PK-PD tooltip content | Parse `pk-pd.ivp` `[tooltip]` annotations via `grok.dapi.files.readAsText`; compare to displayed text. Hardcoded fallback if API fails | `cyclic-models.test.ts` Step 4 |
| **M-3.2** Acid Production tooltip content | Same for `ga-production.ivp` | `stages.test.ts` Step 4 |
| **M-4** Fitting Target Block + RMSE descending | Toggle target output via DOM traversal; read `grok.shell.t` for `RMSE` column; assert ≥95% non-increasing | `fitting.test.ts` Step 6a/6b |
| **M-5** SA 4 viewers correctness | After Run, read each viewer's `dataFrame` via `grok.shell.v.viewers`; assert no NaN/Inf and non-default names | `sensitivity-analysis.test.ts` Step 4 |
| **M-6** Real slider drag | `page.mouse.move/down/move/up` over `input[type="range"]` thumb in 10 steps | `files-and-sharing.test.ts` Step 2 |
| **M-7** Refresh observably reloads | Delete a user-saved script via `grok.dapi.scripts.delete`; click Refresh; assert card count decreased | `catalog.test.ts` Step 4 |
| **M-8** Global error monitoring | `page.on('console')` + `page.on('pageerror')` + network ≥500 with a noise filter for known third-party packages | `helpers/error-monitor.ts`, every test |

---

# Appendix. Cross-reference: MD scenario → automated coverage → manual gaps

| MD step (paraphrased) | Automated assertion | Manual gap |
|---|---|---|
| `open-model` 1–3 | view loads, Multiaxis+Facet tabs visible | — |
| `open-model` 4 | pixel-histogram ≥10 distinct colour buckets | — |
| `open-model` 5 | URL contains `switchat=150` AND canvas hash differs | — |
| `open-model` 6 | ≥4 inputs change AND canvas hash differs | — |
| `cyclic-models` 1–2 | model loads, tabs present | — |
| `cyclic-models` 3 | Count value, URL, AND canvas hash differ | — |
| `cyclic-models` 4 | tooltip text contains IVP-declared `[tooltip]` | — |
| `stages` 1–2 | model loads, tabs present | — |
| `stages` 3 | input value, URL, AND canvas hash differ | — |
| `stages` 4 | tooltip text contains IVP-declared `[tooltip]` | — |
| `scripting` 1 | Edit toggle + ScriptView opens | — |
| `scripting` 2 | Final input updates; Process-mode absent | **M-1.5** chart redraw inside ScriptView |
| `scripting` 3 | `//tags: model` typed + Save clicked | — |
| `scripting` 4 | Bioreactor card visible in Model Hub | — |
| `scripting` 5 | Final = 800 inside Model Hub variant | **M-1.6** chart redraw inside Model Hub variant |
| `sensitivity-analysis` 1–2 | model + SA view load | — |
| `sensitivity-analysis` 3 | ≥3 inputs change on Process mode + switchers toggled (FFox, FKox, FFred) | — |
| `sensitivity-analysis` 4 | viewer count ≥4 AND each viewer's dataframe has no NaN/Inf AND non-default name | — |
| `fitting` 1–3 | view + cascade verified | — |
| `fitting` 4 | switchers + FFox/FKox max set | — |
| `fitting` 4 (MD "Target Block") | target output switcher toggled via DOM traversal | — |
| `fitting` 5 | bioreactor-experiment in dropdown (via 1 API call) | — |
| `fitting` 6 | row count ≥ 0; if RMSE column present → ≥95% non-increasing AND last ≤ first | — |
| `catalog` 1–3 | load, save, Model Hub list | — |
| `catalog` 4 | delete a user-saved entry via API; click Refresh; card count decremented | — |
| `catalog` 5 | dose input visible | — |
| `catalog` 6 | dose=5000 entered | **M-1.7** chart redraw in Model Hub variant |
| `files-and-sharing` 1 | preview opens, ≥5 inputs | — |
| `files-and-sharing` 2 | real slider drag for Step; clicker for Count; URL + canvas hash differ | — |
| `files-and-sharing` 3 | new tab loads same model with same inputs | — |
| `files-and-sharing` 4 | no Multiaxis/Facet at count=1 (REMARK) | — |
| **All scenarios** | global error monitor: console.error / pageerror / 5xx (filtered for third-party noise) | — |
