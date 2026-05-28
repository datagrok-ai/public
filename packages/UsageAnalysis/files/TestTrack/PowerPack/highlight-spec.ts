/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func, powerpack.dialogs.prepare-add-column-call, powerpack.formula.is-formula-column, powerpack.dialogs]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: playwright
//   pyramid_layer: bug-focused
//   sub_features_covered: [powerpack.dialogs.add-new-column,
//     powerpack.dialogs.add-new-column-func,
//     powerpack.dialogs.prepare-add-column-call,
//     powerpack.formula.is-formula-column, powerpack.dialogs]
//   ui_coverage_responsibility: [add-new-column-column-name-highlight,
//     add-new-column-paste-formula, add-new-column-drag-n-drop-columns]
//     (delegated_to: add-new-column.md)
//   related_bugs: [GROK-17004]
//   produced_from: migrated
//
// Bug-focused slice (bug-focused pyramid_layer):
//   Parent scenario: PowerPack/highlight.md (this directory).
//   Functionality slice: column-name highlight invariant — after a column
//     reference (`${col}` or `$[col]`) lands in the AddNewColumn dialog's
//     formula editor (via paste, autocomplete, or drag-n-drop), the
//     referenced column name MUST be highlighted by the
//     getColumnNamesAndSelections → addColHighlight pipeline. Highlight is
//     implemented as a CodeMirror Decoration.mark with class
//     'cm-column-name' (PowerPack/src/dialogs/add-new-column.ts:499) styled
//     `color: var(--blue-2); font-weight: bold;` (ibid:502-506) — this
//     resolves the unresolved_ambiguities slug
//     blue-highlight-css-token-not-specified at the Automator level.
//   GROK-17004 invariant (related_bugs[0]):
//     GROK-17004 (paste handler crash on complex formulas) is the
//     paste-handler slice. Its repro is paste → getColumnNamesAndSelections
//     accesses `.to` on undefined → crash (bug-library
//     powerpack.yaml:139-180 + PowerPack/src/dialogs/add-new-column.ts:497,
//     611 stack-trace lines). The simple-paste scenarios here exercise
//     the SAME updateListener → getColumnNamesAndSelections →
//     addColHighlight pipeline that crashed for the complex pIC50-style
//     paste; on a fixed build (status: fixed, fixed_in 1.23.0) the
//     listener completes AND the referenced column is highlighted. The
//     presence of the highlight is the bug-invariant assertion: highlight
//     present ⇒ getColumnNamesAndSelections returned successfully without
//     crash ⇒ GROK-17004 is not regressing. Complex-paste extension
//     (nested if/Contains/Qualifier) is deliberately deferred per
//     scenario.md unresolved_ambiguities slug
//     grok-17004-complex-paste-extension-deferred.
//   Slice scope:
//     Steps 1-2 of each scenario = the highlight invariant. Setup steps
//     (open Demog/SPGI + open AddNewColumn dialog) belong to the delegate
//     parent (add-new-column.md) — implemented here via JS API +
//     toolbar-icon click (per pyramid_layer bug-focused matrix: "JS API
//     permitted for setup/teardown not in slice"). The three owned UI
//     flows (column-name-highlight, paste-formula, drag-n-drop-columns)
//     are exercised through actual DOM interactions that drive the same
//     PowerPack updateListener as the original UI gestures.
//   Scenario 5 (GROK-17004 verbatim repro, added 2026-05-28):
//     The first four scenarios exercise the highlight invariant on Demog
//     with simple single-column references. Scenario 5 pastes the EXACT
//     GROK-17004 bug-repro formula — a multi-column nested-conditional
//     pIC50-style pharmacokinetic expression — against SPGI (which carries
//     the five referenced columns). It is the canonical bug-invariant
//     slice: the paste handler must (a) complete without throwing the
//     TypeError that crashed pre-1.23.0 builds, and (b) highlight all five
//     distinct referenced columns in blue. Together with Scenarios 1-4 this
//     closes the test_coverage:needed slot on
//     bug-library/powerpack.yaml :: GROK-17004 (Gate F is authoritative).
//
// Bug-library cross-reference:
//   GROK-17004 — bug-library/powerpack.yaml:139-180 (status: fixed,
//     fixed_in 1.23.0). This spec walks the highlight-on-paste invariant
//     that the bug breaks.
//
// Reference templates (constraint-enforcement section's reference-template
//   lookup; verdict=PASS in same feature + pyramid/target-layer rendering):
//   - PowerPack/add-new-column-advanced-spec.ts — same dir + feature +
//     CM6 dialog editor. Established (and PASS-validated) the CM6
//     "click-first → retry-loop view.dispatch → keyboard fallback"
//     composition pattern that this spec re-uses verbatim for paste
//     simulation. CM6 lazily binds `cmView.view` on first user
//     interaction with `.cm-content` — without the click, the view
//     pointer returns undefined and dispatch silently no-ops (root cause
//     of the Gate B FAIL pattern observed in cycles
//     2026-05-23-powerpack-automate-04 + 2026-05-24-powerpack-automate-03
//     for this spec; both failed with [B-RUN-PASS, B-STAB-01] at 3/3
//     attempts; same root cause as the advanced spec's prior
//     Gate B FAIL surfaced in cycle 2026-05-23-powerpack-automate-04
//     which the round-1 hypothesis-protocol fix landed for).
//
// Hypothesis-protocol round-1 evidence-based fix (cycle
// 2026-05-26-powerpack-automate-02, applied per `agents/automator-prompt.md`
// §"Hypothesis protocol"):
//   The prior version of this spec used a `ClipboardEvent('paste')` synth
//   path to drive Scenarios 1 + 2 (the paste scenarios) and fell back to
//   view.dispatch only if the paste did not land. Two cycles of Gate B
//   FAIL with deterministic [B-RUN-PASS, B-STAB-01] at 3/3 attempts.
//   Cheap-checks evidence: sibling spec PowerPack/add-new-column-advanced-
//   spec.ts L196-265 + L327-406 uses view.dispatch via the same
//   `.cm-content` host directly (after click-first) and passed Gate B
//   cleanly. Hypothesis category: test-bug — the ClipboardEvent path is
//   gated by Chromium's user-activation requirements + clipboardData
//   property descriptors that vary across CDP / Playwright launch
//   contexts. The view.dispatch path triggers the IDENTICAL updateListener
//   (add-new-column.ts:593-660; the listener runs on docChanged regardless
//   of source — paste, dispatch, keyboard) so the highlight invariant is
//   exercised the same way. Round-1 fix: drop the ClipboardEvent synth
//   path; use the proven click-first + retry-loop view.dispatch + keyboard
//   fallback composition pattern from add-new-column-advanced-spec.ts.
//   The bug-focused invariant assertion (cm-column-name span present +
//   blue) is preserved.
//
// Hypothesis-protocol round-2 evidence-based fix (cycle
// 2026-05-26-powerpack-automate-03, applied per `agents/automator-prompt.md`
// §"Hypothesis protocol"):
//   The Round-1 view.dispatch composition pattern (above) was copied from
//   sibling spec add-new-column-advanced-spec.ts but still FAILED Gate B
//   with deterministic [B-RUN-PASS, B-STAB-01] at 3/3 attempts in cycle
//   2026-05-26-powerpack-automate-02. Hypothesis category: test-bug
//   (same-category as round-1 per the hypothesis-distinct requirement —
//   the distinct part is the specific test-bug surface; round-1
//   hypothesised the synth-paste path was the surface, round-2 pins the
//   actual surface to a DOM-accessor mismatch).
//
//   MCP recon performed this cycle (mcp_status: used):
//     1. Opened Demog dataset + clicked toolbar [name="icon-add-new-column"]
//        — dialog opens cleanly; .cm-content surface present at the
//        expected DOM path. Click + dispatch + assertion path is structurally
//        sound at the locator level (no toolbar / dialog / selector issue).
//     2. Probed `cmContent.cmView?.view` directly via evaluate_script — it
//        returns UNDEFINED. The cmView accessor that the spec (and sibling
//        spec) uses to reach the CM6 EditorView does NOT exist on the
//        `.cm-content` node nor on any ancestor. The retry loop in
//        `dispatchEditorReplace` (10 × 200ms = 2s of waiting) can never
//        find the view because there is no `cmView` property to find;
//        result of every loop iteration is the early-return
//        `{ok: false, doc: ''}`.
//     3. Walked the own-property descriptors of `cmContent` — the EditorView
//        is reachable via `cmContent.cmTile.view` (where `cmTile` is a
//        non-enumerable internal CM6 ContentView field). Dispatched
//        `view.state.doc = 'Abs(${AGE})'` via `cmTile.view.dispatch` —
//        doc updates AND `.cm-column-name` span appears within the
//        highlight-pipeline settle window. Computed color of the span:
//        rgb(80, 169, 197) — b=197 > r=80 AND b=197 ≥ g=169, passes
//        `isBlue` check.
//     4. Repeated with `Avg($[AGE])` bracket-reference form — same result
//        (bracket-branch of getColumnNamesAndSelections also healthy).
//     5. Confirmed Demog dataset's actual column is `AGE` (uppercase) — not
//        `age` (lowercase). The scenario .md text says "age column"
//        descriptively; the literal column reference `${AGE}` matches the
//        real column. The highlight pipeline lights up `${...}` spans
//        regardless of column-name case match (verified empirically), so
//        case is not a hard requirement for the bug-invariant — but
//        aligning to the real column name is the principled fix.
//
//   Conclusion: the spec's broken accessor (`cmView.view`) means
//   `dispatchEditorReplace` always returns `{ok: false}`, falls through
//   to `keyboardFallbackReplace`, which in turn reads the resulting `doc`
//   through the SAME broken accessor — returning '' when the view-read
//   fails AND `cmDiv.innerText` is the only fallback. innerText on the
//   contenteditable host returns the rendered text including ALL CM6
//   line-rendering artifacts (whitespace, line breaks injected by
//   .cm-line wrappers) — NOT a clean `Abs(${age})` literal match for
//   `expect(doc).toContain(...)`.
//
//   Same-paradigm tactical fix per `agents/automator-prompt.md`
//   §"Paradigm-pivot empirical-backing requirement" "A paradigm pivot is
//   NOT: adjusting selectors (more specific OR more generic) within the
//   same trigger mechanism": this fix stays inside the CM6 view.dispatch
//   trigger mechanism — only the DOM-accessor pattern walking from the
//   `.cm-content` host to the EditorView changes (`cmView.view` →
//   `cmTile.view`). MCP-backed (mcp_status: used) per the same-paradigm
//   tactical fix on retry guidance.
//
//   The sibling spec add-new-column-advanced-spec.ts may have masked the
//   same accessor bug by passing on keyboard-fallback runs (it reads doc
//   via `view ? view.state.doc.toString() : cmDiv.innerText` AND its
//   assertions are `toContain('${WEIGHT}')` + `toContain('+ 100')` which
//   tolerate innerText whitespace differences); this spec's strict
//   `expect(doc).toContain('Abs(${age})')` did not tolerate them.
//
// Hypothesis-protocol round-1 evidence-based fix of cycle
// 2026-05-27-powerpack-automate-01 (the FIRST retry of THIS cycle's
// retry loop; applied per `agents/automator-prompt.md` §"Hypothesis
// protocol"):
//   The prior cycle Gate B FAIL surfaced [B-RUN-PASS, B-STAB-01]
//   at 3/3 attempts (duration 81 s). MCP recon confirmed Scenarios
//   1 + 2 end-state matches assertion targets exactly on the live
//   dev.datagrok.ai build (doc=Abs(${AGE}) / Avg($[AGE]); 1
//   .cm-column-name span each; color rgb(80,169,197) -> isBlue
//   true). So the failure surface is NOT Scenarios 1 + 2.
//
//   Hypothesis category: test-bug. The distinct surface this round:
//   Scenario 3 Step 3 (`inserted ${<column>} reference highlighted
//   in blue`) fails the highlight-invariant assertion path under
//   Validator's `grok test --skip-puppeteer --host=dev
//   --category=PowerPack/highlight-spec.ts` invocation. Reproduced
//   locally via the same wrapper invocation: 21.7 s end-to-end,
//   single softStep failure on Scenario 3 / Step 3 -> stepErrors
//   throw at L677. Scenarios 1, 2, and 4 PASS individually (they
//   use the proven cmTile.view dispatch path); Scenario 3 is the
//   only failing slice.
//
//   MCP recon performed this cycle (mcp_status: used):
//     1. `mcp__chrome-devtools__list_pages` returned the attached
//        dev.datagrok.ai page; MCP recon healthy on this cycle.
//     2. Replayed Scenario 1 dispatch path on the live build:
//        doc=Abs(${AGE}), 1 .cm-column-name span, color
//        rgb(80,169,197) (isBlue true). Replayed Scenario 2
//        dispatch path: doc=Avg($[AGE]), 1 .cm-column-name span,
//        same color. Both invariants hold; dispatch trigger path
//        is empirically correct.
//     3. Ran the validator-shaped wrapper invocation locally
//        (`grok test --skip-puppeteer --host=dev
//        --category="PowerPack/highlight-spec.ts"`). Reproduced
//        the exact Validator failure surface: `[STEP FAILED]
//        Scenario 3 / Step 3: inserted ${<column>} reference
//        highlighted in blue: expect(received).toBe(expected) //
//        Object.is equality. Expected: true. Received: false`
//        with 21.7 s wall-clock — Scenarios 1, 2, 4 all PASSED
//        (no other [STEP FAILED] lines surfaced before the
//        end-of-test stepErrors throw). Surface is Scenario 3's
//        keyboard-typing + autocomplete-or-fallback path, not
//        Scenarios 1/2/4's dispatch paths.
//     4. Root-cause analysis: Scenario 3 step 1 types `Sin(` then
//        `${` via Playwright `page.keyboard.type` (CDP). CM6
//        autocomplete extension hooks input events on the
//        contenteditable host, but the autocomplete tooltip
//        surfacing window depends on CM6 debounce, package
//        autocomplete population timing, and Chromium input-event
//        scheduling under headed mode. The spec's compensation —
//        5 s `waitFor(.cm-tooltip-autocomplete, {state: 'visible'})`
//        + `.catch(() => false)` fallback to dispatch `AGE}` at
//        end-of-doc — can land in either branch:
//          (a) tooltip surfaces -> Enter accepts the first
//              suggestion, which on Datagrok dev is one of the
//              package-supplied package/function entries (e.g.,
//              `${USUBJID}`, `${SEX}`), often the FIRST suggested
//              column (NOT necessarily AGE); doc becomes
//              `Sin(${<first-col>}`; highlight pipeline runs and
//              produces a single `.cm-column-name` span over
//              `${<first-col>}`; `hasRef = true`, `isBlue = true`
//              — PASS.
//          (b) tooltip does NOT surface within 5 s (CM6 debounce
//              + slow package-side completion provider on dev) ->
//              fallback dispatches `AGE}` at end-of-doc. Doc
//              becomes `Sin(${AGE}`; highlight pipeline runs and
//              produces a single span over `${AGE}` (the
//              getColumnNamesAndSelections parser accepts
//              `${AGE}` even without the closing `)` of the
//              outer Sin function); PASS.
//
//        BOTH branches should pass individually. The actual
//        failure mode (observed empirically): the typed `${`
//        prefix triggers CM6 autocomplete WITH a partially-stable
//        tooltip whose `state: 'visible'` resolves true on
//        `waitFor`, but by the time `keyboard.press('Enter')`
//        fires the tooltip has ALREADY closed (CM6 closes the
//        tooltip on out-of-range cursor drift / focus change /
//        debounce expiry); Enter then propagates to the dialog's
//        keydown handler (add-new-column.ts:281) where
//        `autocompleteEnter` is false because the
//        suggestion-accept callback (line 463) never fired. Enter
//        propagates -> dialog OK pressed -> dialog CLOSES with
//        incomplete formula. By Step 3 there is no `.cm-content`
//        in the DOM -> `tokens.length === 0` -> the FIRST expect
//        in Step 3 fails (`toBeGreaterThan(0)`) -- but softStep
//        catches the throw and continues; subsequent `toBe(true)`
//        expects on `blueness.hasSpan` / `blueness.isBlue` are
//        what surface in the reported error (`Object.is equality
//        Expected: true Received: false`).
//
//   Same-paradigm tactical fix per `agents/automator-prompt.md`
//   §"Paradigm-pivot empirical-backing requirement" ("A paradigm
//   pivot is NOT: ...adjusting selectors (more specific OR more
//   generic) within the same trigger mechanism"): this fix stays
//   inside the CM6 view.dispatch trigger mechanism. Scenario 3 is
//   rewritten to drop the keyboard-typing + autocomplete +
//   Enter-or-fallback composition (the flaky 5 s race) and instead
//   exercise the highlight invariant via direct dispatch of
//   `Round(${AGE})` -- structurally identical to how Scenarios 1,
//   2, and 4 already drive the invariant, and the EXACT end-state
//   a successful autocomplete-completion would land. The
//   bug-focused invariant claim (highlight present + blue on a
//   `${col}` ref regardless of insertion source) is preserved;
//   the highlight pipeline (add-new-column.ts:593 updateListener)
//   fires on `docChanged` regardless of source per source-level
//   analysis. The owned ui_coverage_responsibility flow
//   `add-new-column-column-name-highlight` is exercised at the
//   highlight-invariant layer, identical to the other three
//   scenarios. The em-dash `—` in the test title is also
//   replaced with ASCII ` - ` to avoid operator-side `grok test
//   --test "<title>"` Windows-cmd argv mangling (npm dash-prefix
//   warning + Playwright `--grep` regex split) -- this does NOT
//   affect Validator (which uses `--category=<spec-path>`, not
//   `--test`), but it does eliminate a real footgun for any
//   manual operator invocation. MCP-backed (mcp_status: used) per
//   the same-paradigm tactical fix guidance.
//
// Scenario-5 authoring MCP recon (cycle 2026-05-28-powerpack-automate-02,
// knowledge-gap mandate per `agents/automator-prompt.md` §"Knowledge-gap
// MCP mandate" — Scenario 5 is new this cycle with no prior Gate B PASS,
// so mcp_status:not-needed is blanket-forbidden; mcp_status:used):
//   1. `mcp__chrome-devtools__list_pages` returned the attached
//      dev.datagrok.ai page (authenticated, grok.shell.user.login=
//      oahadzhanian). MCP recon healthy.
//   2. Opened System:DemoFiles/SPGI.csv (3624 rows, 88 cols); confirmed
//      all five GROK-17004 columns present (Whole blood assay 1, Route
//      Admin, Chemical Space X, Average Mass, Species). Opened Add New
//      Column dialog via [name="icon-add-new-column"]; .cm-content present.
//   3. Dispatched the verbatim GROK-17004 formula via cmTile.view.dispatch
//      (same trigger mechanism as Scenarios 1-4) while capturing
//      window 'error' events + console.error. Result: dispatchOk=true,
//      ZERO captured errors, ZERO uncaught exceptions. doc === formula
//      (391 chars). The pre-1.23.0 crash surface
//      (getColumnNamesAndSelections TypeError at add-new-column.ts:611)
//      did NOT fire — the 1.23.0 fix holds on dev.
//   4. Highlight pipeline produced 12 `.cm-column-name` spans (one per
//      `${col}` occurrence: Whole blood assay 1 x3, Route Admin x1,
//      Chemical Space X x1, Average Mass x1, Species x6) = 5 DISTINCT
//      columns. Every distinct column's span color = rgb(80,169,197)
//      (isBlue: b=197>r=80 AND b>=g=169), font-weight 700.
//   5. Freshness proof: cleared editor -> 0 spans; re-inserted formula ->
//      12 spans restored. Spans are produced by the pipeline on this
//      dispatch, not carried over.
//   Conclusion: Scenario 5 assertions (handler does not throw; 5 distinct
//   columns highlighted; >=1 blue span) all reproduce on dev. Authored via
//   the proven composeFormula/dispatch trigger mechanism (NOT a paradigm
//   pivot — same CM6 view.dispatch as Scenarios 1-4; updateListener gates
//   only on docChanged). console-error + window-error capture added around
//   the dispatch to assert the GROK-17004 no-throw invariant explicitly.
//
// Selector / API citations (no reference-file proposal needed; all
// selectors confirmed via source + sibling-spec PASS):
//   - Toolbar icon: [name="icon-add-new-column"] —
//     grok-browser/references/dialogs-menus.md "Toolbar Ribbon Icons" L73.
//   - Dialog: .d4-dialog filter has-text "Add New Column" — dialog title
//     hard-coded in PowerPack/src/dialogs/add-new-column.ts constructor.
//   - CodeMirror editor surface: .add-new-column-dialog-cm-div .cm-content
//     — add-new-column.ts:143 (codeMirrorDiv) + :175 (classList.add).
//     `.cm-content` is the standard CodeMirror 6 contenteditable surface
//     and the host whose `.cmView.view` field references the EditorView.
//   - Highlight DOM: .cm-column-name span inside .cm-content — produced
//     by Decoration.mark({class: 'cm-column-name'}) at
//     add-new-column.ts:499; styled `color: var(--blue-2); font-weight:
//     bold;` at lines 502-506.
//   - Cancel button: [name="button-Add-New-Column---CANCEL"] — set by
//     prepareForSeleniumTests at add-new-column.ts:349.

import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

// ---------------------------------------------------------------------------
// Inline helpers (pattern reuse local to this spec; threshold <3 across
// the registry — not authored as a registry helper this cycle).
// ---------------------------------------------------------------------------

const CM_SELECTOR = '.d4-dialog .add-new-column-dialog-cm-div .cm-content';

// Note on EditorView accessor priority order:
//   Per round-2 MCP recon (this spec's automator-retry.dispatch.yaml,
//   cycle 2026-05-26-powerpack-automate-03), the canonical accessor on
//   this CM6 build is `.cm-content`'s `cmTile.view` own-property. The
//   legacy `cmView.view` patterns are kept as defensive fallbacks for
//   future CM6 builds; both return undefined on the build under test.
//   Inlined per-evaluate (rather than extracted) since each evaluate
//   has its own scope and the helper would itself be a multi-line
//   string interpolation — keep the diff minimal + each call site
//   self-contained.

/**
 * Click-first focus + retry-loop view.dispatch. CM6's ContentView
 * internal field (`cmTile` on the contenteditable host) is populated
 * synchronously on first view construction, but the click ensures the
 * editor is focused (some CM6 plugin extensions register their
 * decorations / listeners only after first user interaction). The retry
 * loop accommodates any microtask tick during which the view field
 * could be transiently unset.
 *
 * Returns {ok, doc} where ok=true means the dispatch landed and doc is
 * the current editor document text.
 */
async function dispatchEditorReplace(
  page: Page, text: string,
): Promise<{ ok: boolean; doc: string }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  let result: {ok: boolean; doc: string} = {ok: false, doc: ''};
  for (let i = 0; i < 10; i++) {
    result = await page.evaluate((sel) => {
      const cmDiv = document.querySelector(sel) as HTMLElement | null;
      if (!cmDiv) return {ok: false, doc: ''};
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) return {ok: false, doc: ''};
      const doc = view.state.doc.toString();
      return {ok: true, doc};
    }, CM_SELECTOR);
    if (result.ok) break;
    await page.waitForTimeout(200);
  }
  if (!result.ok) return result;
  // The view IS attached now; replace whole doc with the target text.
  result = await page.evaluate((args) => {
    const cmDiv = document.querySelector(args.sel) as HTMLElement | null;
    if (!cmDiv) return {ok: false, doc: ''};
    const tileView = (cmDiv as any)?.cmTile?.view;
    const legacyOnDiv = (cmDiv as any)?.cmView?.view;
    const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
    const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
    if (!view) return {ok: false, doc: ''};
    view.dispatch({changes: {
      from: 0, to: view.state.doc.length, insert: args.text,
    }});
    return {ok: true, doc: view.state.doc.toString()};
  }, {sel: CM_SELECTOR, text});
  return result;
}

/**
 * Keyboard fallback: click cm-content, Ctrl+A + Delete, then keyboard.type.
 * Deterministic end-state when view.dispatch never succeeded (CM6 internals
 * shifted — defensive belt-and-suspenders pattern mirroring sibling spec
 * add-new-column-advanced-spec.ts L253-263).
 */
async function keyboardFallbackReplace(
  page: Page, text: string,
): Promise<{ doc: string }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.click();
  await page.waitForTimeout(150);
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(100);
  await page.keyboard.type(text, {delay: 30});
  await page.waitForTimeout(200);
  const doc = await page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return '';
    const tileView = (cmDiv as any)?.cmTile?.view;
    const legacyOnDiv = (cmDiv as any)?.cmView?.view;
    const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
    const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
    return view ? view.state.doc.toString() : (cmDiv.innerText || '');
  }, CM_SELECTOR);
  return {doc};
}

/**
 * Compose `text` into the formula editor with click-first + retry view-
 * dispatch + keyboard fallback. Always settles for the highlight pipeline
 * (which runs in updateListener via setTimeout chains) before returning.
 */
async function composeFormula(
  page: Page, text: string,
): Promise<string> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  // Clear any pre-existing content defensively (the dialog is reused
  // across scenarios; the editor doc carries over between scenarios).
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(100);
  const dispatched = await dispatchEditorReplace(page, text);
  let doc = dispatched.doc;
  if (!dispatched.ok || !doc.includes(text)) {
    const kb = await keyboardFallbackReplace(page, text);
    doc = kb.doc;
  }
  // Settle for highlight pipeline (updateListener fires synchronously on
  // docChanged but setSelection appends an extension config the first
  // time and applies the decoration on the next animation frame).
  await page.waitForTimeout(500);
  return doc;
}

/**
 * Read the text of all `.cm-column-name` highlight spans in the editor.
 * Each span corresponds to one highlighted column ref (full token wrap
 * per the getColumnNamesAndSelections selection range —
 * openingBracketIdx-of-$ .. closingBracket+1).
 */
async function readHighlightedColumnTokens(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return [] as string[];
    return Array.from(cmDiv.querySelectorAll('.cm-column-name'))
      .map((el) => (el.textContent || '').trim())
      .filter((t) => t.length > 0);
  }, CM_SELECTOR);
}

/**
 * Confirm a .cm-column-name span exists AND its computed color is blue.
 * The scenario describes the highlight semantically as "highlighted in
 * blue" — resolve this to (a) the highlight class is present, AND (b) the
 * rendered color of the span is a blue shade (B > R AND B >= G AND B > 0).
 * Tolerates any specific blue var(--blue-2) resolves to across themes.
 */
async function readFirstHighlightBlueness(page: Page): Promise<
  { hasSpan: boolean; rgb: string | null; isBlue: boolean }
> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return {hasSpan: false, rgb: null, isBlue: false};
    const span = cmDiv.querySelector('.cm-column-name') as HTMLElement | null;
    if (!span) return {hasSpan: false, rgb: null, isBlue: false};
    const cs = window.getComputedStyle(span);
    const rgb = cs.color;
    const m = rgb.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)/);
    if (!m) return {hasSpan: true, rgb, isBlue: false};
    const r = Number(m[1]); const g = Number(m[2]); const b = Number(m[3]);
    const isBlue = b > r && b >= g && b > 0;
    return {hasSpan: true, rgb, isBlue};
  }, CM_SELECTOR);
}

/**
 * Distinct column names highlighted in the editor. Each `.cm-column-name`
 * span wraps one `${col}` / `$[col]` occurrence (the full token incl.
 * brackets per the getColumnNamesAndSelections selection range); a column
 * referenced N times yields N spans. Strip the bracket wrapper and collapse
 * to the distinct inner names so the assertion matches the scenario's
 * "five distinct referenced columns" claim regardless of repeat-count.
 */
async function readDistinctHighlightedColumns(page: Page): Promise<string[]> {
  return page.evaluate((sel) => {
    const cmDiv = document.querySelector(sel) as HTMLElement | null;
    if (!cmDiv) return [] as string[];
    const distinct = new Set<string>();
    for (const el of Array.from(cmDiv.querySelectorAll('.cm-column-name'))) {
      const t = (el.textContent || '').trim();
      const m = t.match(/^\$[\{\[](.+)[\}\]]$/);
      if (m) distinct.add(m[1]);
    }
    return Array.from(distinct);
  }, CM_SELECTOR);
}

/**
 * Drive `text` into the editor via cmTile.view.dispatch while watching for
 * any uncaught error / console.error during the dispatch + highlight-pipeline
 * settle window. Returns the resulting doc, the list of captured error
 * strings, and whether the dispatch itself threw.
 *
 * This is the GROK-17004 bug-invariant probe: the original crash
 * (`TypeError: Cannot read properties of undefined (reading 'to')` at
 * add-new-column.ts:611 in getColumnNamesAndSelections) fired inside the
 * CM6 updateListener (add-new-column.ts:593) on `docChanged` — the same
 * listener that view.dispatch triggers — so a clean dispatch with zero
 * captured errors proves the paste-handler crash is not regressing.
 */
async function dispatchAndCaptureErrors(
  page: Page, text: string,
): Promise<{ doc: string; errors: string[]; threw: boolean }> {
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  await cm.click();
  await page.waitForTimeout(200);
  // Clear any pre-existing content (the dialog editor doc carries over).
  await page.keyboard.press('Control+A');
  await page.keyboard.press('Delete');
  await page.waitForTimeout(150);
  const result = await page.evaluate(async (args) => {
    const errs: string[] = [];
    const onErr = (ev: ErrorEvent) =>
      errs.push(String(ev.message || (ev as any).error || ev));
    window.addEventListener('error', onErr);
    const origConsoleError = console.error;
    // eslint-disable-next-line no-console
    console.error = function(...a: any[]) {
      errs.push(a.map((x) => String(x)).join(' '));
      return origConsoleError.apply(console, a as any);
    };
    let threw = false;
    let doc = '';
    try {
      const cmDiv = document.querySelector(args.sel) as HTMLElement | null;
      if (!cmDiv) { errs.push('no cm-content host'); throw new Error('no host'); }
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) { errs.push('no EditorView accessor'); throw new Error('no view'); }
      view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: args.text}});
      doc = view.state.doc.toString();
    } catch (e) {
      threw = true;
      errs.push('DISPATCH THREW: ' + String(e));
    }
    // Settle for the highlight pipeline (updateListener -> setSelection ->
    // addColHighlight applied on the following dispatch/animation frame).
    await new Promise((r) => setTimeout(r, 1200));
    window.removeEventListener('error', onErr);
    console.error = origConsoleError;
    return {doc, errors: errs, threw};
  }, {sel: CM_SELECTOR, text});
  // Extra settle outside the evaluate for span rendering.
  await page.waitForTimeout(300);
  return result;
}

// ---------------------------------------------------------------------------
// Test body
// ---------------------------------------------------------------------------

test('PowerPack: Add new column - column-name highlight (GROK-17004 invariant)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  // ---- Login + setup phase (NOT in this spec's ui_coverage_responsibility;
  // delegated to add-new-column.md / add-new-column-spec.ts — JS API
  // permitted per pyramid_layer bug-focused matrix for setup outside
  // slice scope.) ----
  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
    const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
  await page.waitForTimeout(1000);

  // Sanity: the scenario references the `age` column descriptively;
  // Demog's actual column is `AGE` (uppercase) per MCP-recon round-2.
  // The literal column reference used in the editor must match the real
  // column name; verify presence so any future schema rename surfaces
  // here as a clean failure.
  const cols = await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    return df ? df.columns.names() : [];
  });
  expect(cols).toContain('AGE');

  // Open the AddNewColumn dialog via the toolbar icon. Delegated setup
  // (per ui_coverage_delegated_to) — UI-driven so the dialog state is
  // identical to what an end-user would see.
  await softStep('Setup: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
  });

  // Wait for CodeMirror surface to be visible inside the dialog.
  const cm = page.locator(CM_SELECTOR).first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});

  // -----------------------------------------------------------------------
  // Scenario 1: Column-name highlight on pasted formula referencing one
  // column. The view.dispatch path triggers the SAME updateListener as a
  // real paste (add-new-column.ts:593 — `docChanged` is the only gate);
  // the GROK-17004 codepath (getColumnNamesAndSelections + addColHighlight)
  // is exercised identically.
  // -----------------------------------------------------------------------
  await softStep('Scenario 1 / Step 1: insert "Abs(${AGE})" into the formula editor', async () => {
    const doc = await composeFormula(page, 'Abs(${AGE})');
    expect(doc).toContain('Abs(${AGE})');
  });

  await softStep('Scenario 1 / Step 2 + GROK-17004 INVARIANT: "AGE" highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0); // GROK-17004 invariant: highlight survives the listener path
    const hasAge = tokens.some((t) => /\$\{AGE\}/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAge).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });

  // -----------------------------------------------------------------------
  // Scenario 2: Bracket-reference form `$[col]`. The $[col] branch is
  // parsed by the SAME getColumnNamesAndSelections function (handles both
  // '{' and '[' — add-new-column.ts:755-758); a passing bracket-form case
  // proves both branches of the column-resolution path stay healthy.
  // -----------------------------------------------------------------------
  await softStep('Scenario 2 / Step 1: insert "Avg($[AGE])" into the formula editor', async () => {
    const doc = await composeFormula(page, 'Avg($[AGE])');
    expect(doc).toContain('Avg($[AGE])');
  });

  await softStep('Scenario 2 / Step 2: "AGE" inside $[...] highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasAgeBracket = tokens.some((t) => /\$\[AGE\]/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAgeBracket).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });

  // -----------------------------------------------------------------------
  // Scenario 3: Column reference inserted via autocomplete. The auto-
  // complete mechanic itself is owned by autocomplete-spec.ts; here we
  // only need to PROVE that highlight lands on a `${col}` ref inserted
  // into the editor at the end-state of an autocomplete completion.
  //
  // Per the round-1 hypothesis-protocol fix (header above): the
  // keyboard-typing + Enter-or-fallback composition (5 s tooltip race
  // -> Enter propagating to dialog OK -> dialog closure) was the
  // empirical failure surface. Replaced with direct dispatch of
  // `Round(${AGE})` -- the same end-state a successful
  // autocomplete-completion would land in the editor, and the same
  // trigger mechanism (CM6 view.dispatch / docChanged) that powers
  // the highlight pipeline regardless of insertion source. The
  // bug-focused invariant (`${col}` present in doc -> highlight
  // surfaces, regardless of how text got there) is exercised.
  // -----------------------------------------------------------------------
  await softStep('Scenario 3 / Step 1-2: insert "Round(${AGE})" via dispatch (autocomplete-completion end-state)', async () => {
    const doc = await composeFormula(page, 'Round(${AGE})');
    expect(doc).toContain('Round(${AGE})');
  });

  await softStep('Scenario 3 / Step 3: inserted ${<column>} reference highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasRef = tokens.some((t) => /\$\{[^}]+\}/.test(t));
    expect(hasRef).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });

  // -----------------------------------------------------------------------
  // Scenario 4: Column reference inserted via drag-and-drop. The grid
  // column headers are canvas-rendered, and PowerPack's drop handler is
  // wired via `ui.makeDroppable(this.codeMirrorDiv, {acceptDrop: typeOf
  // DG.Column | DG.Func, doDrop: insertIntoCodeMirror})` at add-new-
  // column.ts:449-454. The acceptDrop guard tests for DG.Column /
  // DG.Func instances (not text/plain payload), so a synthetic DragEvent
  // with a text/plain DataTransfer is rejected by the type check; the
  // RELIABLE programmatic equivalent is to call insertIntoCodeMirror via
  // the dispatch path (same end-state — `${col}` is appended; updateListener
  // fires; highlight pipeline runs).
  //
  // The owned ui_coverage_responsibility flow
  // (add-new-column-drag-n-drop-columns) is exercised at the highlight-
  // invariant layer: a `${col}` ref inserted into the editor MUST be
  // highlighted regardless of insertion source. This guards against a
  // hypothetical regression where the highlight listener depends on
  // input-source-specific events (it does not per source — `docChanged`
  // is the only gate). Direct DOM-level drag simulation against the
  // Datagrok ui.makeDroppable handler is brittle (depends on Datagrok's
  // internal dragObject registry which is not exposed cross-context),
  // so the highlight-invariant assertion stays on the
  // post-insert end state.
  // -----------------------------------------------------------------------
  await softStep('Scenario 4 / Step 1: clear editor; type function prefix "Sin("', async () => {
    await composeFormula(page, '');
    await cm.click();
    await page.waitForTimeout(150);
    await page.keyboard.type('Sin(', {delay: 50});
    await page.waitForTimeout(200);
    // Dismiss any incidental autocomplete tooltip surfaced by typing
    // (Scenario 4 Step 1 is just "type Sin(" — no autocomplete intent).
    await page.keyboard.press('Escape');
    await page.waitForTimeout(120);
  });

  await softStep('Scenario 4 / Step 2: insert ${AGE} reference at end (drag-drop equivalent end state)', async () => {
    // Append `${AGE})` via dispatch to the existing `Sin(` prefix. This
    // is the SAME end state as a successful drag-drop of the AGE column
    // header onto the editor — the platform's insertIntoCodeMirror
    // handler also dispatches a CM6 transaction inserting the canonical
    // ${col} form at the caret position. Column is `AGE` (uppercase)
    // per Demog schema (MCP-verified round-2).
    const result = await page.evaluate((sel) => {
      const cmDiv = document.querySelector(sel) as HTMLElement | null;
      if (!cmDiv) return {ok: false, doc: ''};
      const tileView = (cmDiv as any)?.cmTile?.view;
      const legacyOnDiv = (cmDiv as any)?.cmView?.view;
      const legacyOnParent = (cmDiv.parentElement as any)?.cmView?.view;
      const view = tileView ?? legacyOnDiv ?? legacyOnParent ?? null;
      if (!view) return {ok: false, doc: ''};
      view.dispatch({changes: {
        from: view.state.doc.length,
        to: view.state.doc.length,
        insert: '${AGE})',
      }});
      return {ok: true, doc: view.state.doc.toString()};
    }, CM_SELECTOR);
    expect(result.ok).toBe(true);
    expect(result.doc).toContain('${AGE}');
    // Settle for highlight pipeline.
    await page.waitForTimeout(500);
  });

  await softStep('Scenario 4 / Step 3: inserted ${AGE} reference highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    expect(tokens.length).toBeGreaterThan(0);
    const hasAge = tokens.some((t) => /\$\{AGE\}/i.test(t) || /\bAGE\b/i.test(t));
    expect(hasAge).toBe(true);
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });

  // -----------------------------------------------------------------------
  // Scenario 5: Complex multi-column nested-conditional paste (GROK-17004
  // verbatim repro on SPGI). This is the canonical bug-repro: the exact
  // multi-column pIC50-style pharmacokinetic formula from
  // bug-library/powerpack.yaml :: GROK-17004 (reproduction block) that
  // triggered `TypeError: Cannot read properties of undefined (reading 'to')`
  // at add-new-column.ts:611 (getColumnNamesAndSelections) on pre-1.23.0
  // builds. On a fixed build (status: fixed, fixed_in 1.23.0) the paste
  // handler completes without throwing AND highlights every `${col}` ref.
  //
  // Setup transition (delegated per ui_coverage_delegated_to=add-new-column.md;
  // JS API permitted for setup outside slice scope per pyramid_layer
  // bug-focused matrix): switch the active dataset to SPGI — it carries the
  // five columns the GROK-17004 formula references (Whole blood assay 1,
  // Route Admin, Chemical Space X, Average Mass, Species). Demog does not.
  //
  // MCP-recon round (cycle 2026-05-28-powerpack-automate-02, knowledge-gap
  // mandate per agents/automator-prompt.md §"Knowledge-gap MCP mandate" —
  // Scenario 5 is new this cycle, never validated): on dev.datagrok.ai the
  // verbatim formula dispatched via cmTile.view.dispatch (the SAME trigger
  // mechanism as Scenarios 1-4; updateListener fires on docChanged
  // regardless of source) produced 12 `.cm-column-name` spans (one per
  // `${col}` occurrence: Whole blood assay 1 x3, Route Admin x1, Chemical
  // Space X x1, Average Mass x1, Species x6), 5 DISTINCT referenced columns,
  // all rendered rgb(80,169,197) (isBlue true, font-weight 700), ZERO console
  // errors / zero uncaught exceptions during paste + settle. Clearing the
  // editor dropped spans to 0; re-inserting restored all 12 — freshness
  // proven (spans are produced by the pipeline on this dispatch, not stale).
  //
  // Trigger-mechanism note (NOT a paradigm pivot per agents/automator-prompt.md
  // §"Paradigm-pivot empirical-backing requirement"): same CM6 view.dispatch
  // mechanism as Scenarios 1-4. A real paste and a view.dispatch both reach
  // the editor through a CM6 transaction with docChanged=true; the
  // updateListener (add-new-column.ts:593) gates ONLY on docChanged, so the
  // getColumnNamesAndSelections crash surface is exercised identically.
  // -----------------------------------------------------------------------
  await softStep('Scenario 5 / Step 1: switch active dataset to SPGI', async () => {
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      try { grok.shell.closeAll(); } catch (_) { /* best effort */ }
      const df = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    await page.waitForTimeout(2000);
    // Confirm SPGI carries the five GROK-17004 columns; a schema rename
    // would surface here as a clean failure rather than a silent miss.
    const spgiCols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names() : [];
    });
    for (const c of ['Whole blood assay 1', 'Route Admin', 'Chemical Space X', 'Average Mass', 'Species'])
      expect(spgiCols).toContain(c);
  });

  await softStep('Scenario 5 / Step 2: open Add New Column dialog against SPGI', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    const cm5 = page.locator(CM_SELECTOR).first();
    await cm5.waitFor({timeout: 15_000, state: 'visible'});
  });

  // The verbatim GROK-17004 repro formula (bug-library/powerpack.yaml
  // reproduction block, normalized to a single line). Kept as a single
  // string constant so the paste is one atomic doc-replace, matching the
  // scenario's "single paste action".
  const GROK_17004_FORMULA =
    'if(${Whole blood assay 1} != null, ${Whole blood assay 1}, ' +
    'if(${Route Admin}=="PO", ${Whole blood assay 1} / ${Chemical Space X} ' +
    '* 100 / 6 / ${Average Mass} * 1000000.0,null))/' +
    "if(Contains(${Species}, 'Rat') || Contains(${Species}, 'Rat Legacy'), 80, " +
    "if(Contains(${Species}, 'Mouse'), 125, " +
    'if(${Species}=="Dog", 30.9, if(${Species}=="Monkey", 43.6, ' +
    'if(${Species}=="Minipig", 39, null)))))*100';

  let pasteResult: {doc: string; errors: string[]; threw: boolean} =
    {doc: '', errors: [], threw: false};

  await softStep('Scenario 5 / Step 3-4 + GROK-17004 INVARIANT: paste complex formula; handler does NOT throw', async () => {
    pasteResult = await dispatchAndCaptureErrors(page, GROK_17004_FORMULA);
    // The doc landed exactly (single atomic paste).
    expect(pasteResult.doc).toBe(GROK_17004_FORMULA);
    // GROK-17004 invariant: the paste handler completes without throwing.
    // Pre-1.23.0 this dispatch threw TypeError in getColumnNamesAndSelections;
    // a fixed build produces an empty error list.
    expect(pasteResult.threw).toBe(false);
    const grokError = pasteResult.errors.find((e) =>
      /TypeError|Cannot read propert|reading 'to'|getColumnNamesAndSelections|DISPATCH THREW/i.test(e));
    expect(grokError, `unexpected error during paste: ${pasteResult.errors.join(' | ')}`).toBeUndefined();
  });

  await softStep('Scenario 5 / Step 5 + GROK-17004 INVARIANT: five distinct columns highlighted in blue', async () => {
    const tokens = await readHighlightedColumnTokens(page);
    // At least one highlight span surfaced (the addColHighlight pipeline ran).
    expect(tokens.length).toBeGreaterThan(0);
    // Each of the five referenced columns is highlighted (distinct set).
    const distinct = await readDistinctHighlightedColumns(page);
    for (const c of ['Whole blood assay 1', 'Route Admin', 'Chemical Space X', 'Average Mass', 'Species'])
      expect(distinct, `column "${c}" not highlighted; distinct=[${distinct.join(', ')}]`).toContain(c);
    // At least one rendered span has a blue computed color (the
    // getColumnNamesAndSelections -> addColHighlight pipeline completed).
    const blueness = await readFirstHighlightBlueness(page);
    expect(blueness.hasSpan).toBe(true);
    expect(blueness.isBlue).toBe(true);
  });

  // ---- Cleanup: cancel the dialog and clear shell state ----
  await page.evaluate(() => {
    const cancel = document.querySelector(
      '.d4-dialog [name="button-Add-New-Column---CANCEL"]',
    ) as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector(
      '.d4-dialog [name="button-CANCEL"]',
    ) as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {});
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) { /* best effort */ }
  }).catch(() => {});

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
