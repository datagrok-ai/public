# Diff Studio — Export Dialog Implementation Plan

## Overview

Add two actions to Datagrok's right-hand **Context panel**, surfaced **only while the Diff Studio app is running with a loaded model**: **Export as Markdown…** and **Export as LaTeX…**. Each opens a modal dialog with a live source preview and a small set of options. The dialog uses Datagrok UI primitives (`ui.dialog`, `ui.input.*`, `ui.splitH`) and the already-used `diff-grok` library for all conversion logic. No `package.ts` changes; no global platform registration.

## Scope decisions (locked in)

- `diff-grok` is already a dependency of the DiffStudio package; import directly as `from 'diff-grok'`.
- Only `convertIvpToLatex` and `ConvertOptions` are used from `diff-grok`. `wrapTexDocument` is **not** part of the library's public top-level entry, so the LaTeX preamble is reproduced locally as a small `wrapTex()` helper inside `export-dialog.ts`.
- **No `ObjectHandler` and no `package.ts` changes.** Actions live on the Context panel only while Diff Studio is active with a model. **Invariant:** whenever Diff Studio is the active view AND a model is loaded, the Context panel always contains a **Diff Studio Export** pane, regardless of what object is currently selected (model, graph, column). The pane is attached at runtime from the `DiffStudio` app class via `grok.events.onAccordionConstructed`, scoped to the app's view lifetime via `view.subs.push(...)`.
- **Standalone document** is a user-controlled toggle (LaTeX only). When on, the preview and the downloaded file both include `\documentclass{article}` + `\usepackage{amsmath, amssymb, booktabs}` + `\begin{document} … \end{document}`. When off, both show the body fragment only. Preview = Download at all times — no silent wrapping on Download.
- Markdown output is always self-contained; the Standalone toggle is hidden / disabled when format is Markdown.
- The ribbon **Download** icon (top of Diff Studio) is **not** touched — it keeps downloading the raw `.ivp`.
- All four **Include** checkboxes default to `on`. The user un-ticks what they do not need. **Standalone document** also defaults to `on`, so first-time LaTeX users get a file that compiles with `pdflatex` out of the box.
- Settings are **not** persisted between sessions. The dialog opens with the same defaults every time.
- Every checkbox, choice input, and button has a tooltip. Tooltip text is taken from the `diff-grok` public documentation (README and `ConvertOptions` JSDoc in `types.ts`) to stay consistent with the library's own vocabulary.

## Stage 1. Wiring `diff-grok`

`diff-grok` is already available in the DiffStudio package. In the dialog module:

```typescript
import {convertIvpToLatex, ConvertOptions} from 'diff-grok';
```

No `package.json` or webpack changes needed.

**Preamble for standalone `.tex`.** `wrapTexDocument` exists in `diff-grok` but is not re-exported from the top-level entry, so we replicate its logic locally (ten lines, no deep imports):

```typescript
// mirrors diff-grok/src/latex-export/output/tex-formatter.ts
// keep packages in sync with what buildLatexDocument actually emits:
// align (amsmath), symbols (amssymb), booktabs rules (booktabs)
function wrapTex(content: string): string {
  return [
    '\\documentclass{article}',
    '\\usepackage{amsmath}',
    '\\usepackage{amssymb}',
    '\\usepackage{booktabs}',
    '',
    '\\begin{document}',
    '',
    content,
    '',
    '\\end{document}',
  ].join('\n');
}
```

If a future `diff-grok` release changes the required packages (e.g. adds `siunitx` or `tikz`), re-check this list when bumping the dependency version.

## Stage 2. Dialog module

Create `src/export/export-dialog.ts` with a single public function:

```typescript
export function showExportDialog(
  ivpText: string,
  modelName: string,
  initialFormat: 'markdown' | 'latex',
): void;
```

Built with `ui.dialog(...)`. Layout: `ui.splitH([optionsPanel, previewPanel])` with a fixed-width left column (~200 px) and a minimum dialog height (~500 px) set via `dialog.root.style`. Cancel button is added by Datagrok automatically; **Copy** and **Download** are explicit `addButton` entries.

The dialog title follows the format: "Export as Markdown" / "Export as LaTeX". It updates when the user switches the **Format** input.

## Stage 3. Options via `ui.input.*`

Every option is a native Datagrok input, so the dialog inherits platform styling and behavior automatically. Tooltips are attached via the standard `tooltipText` field on each input.

The dialog maintains a small state object. `ConvertOptions` covers everything `diff-grok` understands; **Standalone document** is a dialog-local flag — not part of `ConvertOptions` — applied as a post-processing step.

```typescript
const state: Partial<ConvertOptions> & {standalone: boolean} = {
  format: initialFormat,
  includeMetadata: true,
  includeInits: true,
  includeParameters: true,
  includeConstants: true,
  compact: true,      // auto-disabled below for models with >3 equations
  useCdot: true,
  standalone: true,   // LaTeX only; ignored when format === 'markdown'
};

const formatInput = ui.input.choice('Format', {
  value: state.format,
  items: ['markdown', 'latex'],
  tooltipText: 'Output format: LaTeX (.tex) or Markdown with LaTeX math (.md).',
});

const metaInput = ui.input.bool('Title & description', {
  value: state.includeMetadata!,
  tooltipText: 'Include metadata sections: model name, description, and comment.',
});

const initsInput = ui.input.bool('Initial conditions', {
  value: state.includeInits!,
  tooltipText: 'Include the initial conditions table (variable values at the start of the simulation).',
});

const paramsInput = ui.input.bool('Parameters', {
  value: state.includeParameters!,
  tooltipText: 'Include the parameters table (values that generate UI controls for model exploration).',
});

const constsInput = ui.input.bool('Constants', {
  value: state.includeConstants!,
  tooltipText: 'Include the constants table (fixed values used in equations).',
});

const compactInput = ui.input.bool('Compact layout', {
  value: state.compact!,
  tooltipText: 'Compact mode: no section headings, inline initial conditions. Suitable for small models (up to three equations).',
});

const cdotInput = ui.input.choice('Multiplication', {
  value: 'a · b',
  items: ['a · b', 'ab'],
  tooltipText: 'Symbol used for multiplication. "a · b" renders \\cdot; "ab" uses juxtaposition, common in physics papers.',
});

const standaloneInput = ui.input.bool('Standalone document', {
  value: state.standalone,
  tooltipText: 'Wrap the output with \\documentclass{article} and the required \\usepackage lines so the file compiles with pdflatex out of the box. Turn off to get a body fragment for pasting into an existing LaTeX document. Ignored for Markdown.',
});
```

Groupings in the left column: a small `h4`-style label per group (`Format`, `Include`, `Style`) followed by the inputs stacked vertically. **Standalone document** sits in the `Style` group, directly under **Multiplication**. Its input `root` is hidden (`display: none`) whenever `formatInput.value === 'markdown'` — reshown on switch back to LaTeX. No accordion — three short groups fit on screen without collapsing.

Tooltip strings are consciously drawn from the `diff-grok` README (the "What gets converted" and "Compact mode" sections) so a user who reads both the tooltip and the library docs sees consistent wording.

## Stage 4. Live preview

The preview is a monospace `<pre>` element. Regeneration is synchronous — `convertIvpToLatex` is fast enough that neither debouncing nor a spinner is warranted.

```typescript
const previewEl = ui.element('pre');
previewEl.style.cssText =
  'font-family: var(--grok-font-family-monospace, monospace); ' +
  'font-size: 12px; margin: 0; padding: 12px 16px; ' +
  'overflow: auto; white-space: pre-wrap; max-height: 360px;';

let lastResult: string | null = null;

function regenerate(): void {
  try {
    const opts = collectOptions();
    const body = convertIvpToLatex(ivpText, opts);
    const result = (opts.format === 'latex' && standaloneInput.value)
      ? wrapTex(body)
      : body;
    previewEl.textContent = result;
    lastResult = result;
    setActionsEnabled(true);
  } catch (e) {
    previewEl.textContent = `// Error: ${(e as Error).message}`;
    lastResult = null;
    setActionsEnabled(false);
  }
}

[formatInput, metaInput, initsInput, paramsInput,
 constsInput, compactInput, cdotInput, standaloneInput]
  .forEach((input) => input.onChanged.subscribe(regenerate));

// show/hide Standalone based on format
formatInput.onChanged.subscribe(() => {
  standaloneInput.root.style.display =
    formatInput.value === 'latex' ? '' : 'none';
});

regenerate(); // initial render
```

Because `regenerate()` applies `wrapTex` itself, **the preview always shows exactly what Download would save.** No silent wrapping at Download time.

`collectOptions()` maps the seven `ConvertOptions` inputs back into a `Partial<ConvertOptions>`. `standalone` is intentionally excluded — it is consumed locally by `regenerate()`, not by `diff-grok`:

```typescript
function collectOptions(): Partial<ConvertOptions> {
  return {
    format: formatInput.value as 'markdown' | 'latex',
    includeMetadata: metaInput.value,
    includeInits: initsInput.value,
    includeParameters: paramsInput.value,
    includeConstants: constsInput.value,
    compact: compactInput.value,
    useCdot: cdotInput.value === 'a · b',
  };
}
```

## Stage 5. Copy and Download

File name is a separate `ui.input.string` positioned below the preview. It auto-updates its extension when the format changes.

```typescript
const fileNameInput = ui.input.string('File name', {
  value: `${modelName}-${isoDate()}.${initialFormat === 'markdown' ? 'md' : 'tex'}`,
  tooltipText: 'Name of the file created by Download. The extension updates automatically with the selected format.',
});

formatInput.onChanged.subscribe(() => {
  const ext = formatInput.value === 'markdown' ? 'md' : 'tex';
  fileNameInput.value = fileNameInput.value.replace(/\.(md|tex)$/, `.${ext}`);
});
```

Buttons — tooltips set via `ui.tooltip.bind(buttonElement, 'text')`:

```typescript
const dialog = ui.dialog('Export as ' + labelForFormat(initialFormat))
  .add(ui.splitH([optionsPanel, previewPanel]))
  .add(fileNameInput.root)
  .addButton('Copy', async () => {
    if (!lastResult) return;
    await navigator.clipboard.writeText(lastResult);
    grok.shell.info('Copied to clipboard');
    // dialog stays open on Copy
  })
  .addButton('Download', () => {
    if (!lastResult) return;
    // lastResult already reflects the Standalone toggle — no extra wrapping here.
    DG.Utils.download(fileNameInput.value, lastResult);
    dialog.close();
  });

ui.tooltip.bind(copyBtn, 'Copy the generated source to the clipboard, exactly as shown in the preview. The dialog stays open so you can tweak options and copy again.');
ui.tooltip.bind(downloadBtn, 'Save the preview as a file. Use "Standalone document" on the left to control whether the LaTeX output is a compilable file or a body fragment.');

dialog.showModal();
```

If Datagrok's default OK button cannot be suppressed cleanly, use three explicit `addButton` calls (Cancel, Download, Copy) and leave `onOK` empty.

On Copy, the dialog stays open — this matches the common scientific workflow of copying, testing the paste, and coming back to adjust options. On Download, the dialog closes.

## Stage 6. Context-panel integration (app-scoped, no `package.ts` changes)

The two actions live on Datagrok's right-hand **Context panel**. **Invariant:** while Diff Studio is the currently active view *and* a model is loaded, the Context panel always contains a **Diff Studio Export** pane — regardless of what object the user has selected (the model itself, a graph, a column inside an output table). Switch to another view → pane disappears from any freshly built accordion. Close Diff Studio → subscription is torn down.

No `ObjectHandler`, no `//tags: init`, no global side-effects.

### 6.1 Current model, held by the app

`DiffStudio` maintains a single `currentModelRef: DiffStudioModelRef | null` member — rebuilt whenever the loaded model changes (open file, switch tab, apply template, save edit). The ref is a plain value object:

```typescript
// src/export/diff-studio-model-ref.ts
export class DiffStudioModelRef {
  constructor(
    public readonly ivpText: string,
    public readonly modelName: string,
  ) {}
}
```

When no model is loaded (empty editor, initial state before a template/file is picked), `currentModelRef` stays `null` and the pane does not appear. This is the **only** thing that must be kept in sync — the rest is event-driven.

Optionally, the app may also assign `grok.shell.o = this.currentModelRef` at model-load time so the Context panel has something meaningful to show immediately on app open; this is cosmetic and does not affect the pane's presence (the pane attaches independently of `grok.shell.o`).

### 6.2 Subscription: every accordion, app-scoped guard

Inside app startup (`runSolverApp` / `runSolverDemoApp`, wherever `this.view` becomes available), subscribe to `grok.events.onAccordionConstructed` and push the subscription onto `view.subs` for automatic `detach()`-time unsubscribe:

```typescript
// inside DiffStudio app startup, after this.view is created
const EXPORT_PANE_NAME = 'Diff Studio Export';

this.view.subs.push(
  grok.events.onAccordionConstructed.subscribe((acc: DG.Accordion) => {
    // App-scoped guard: only while Diff Studio is active AND has a model loaded.
    if (grok.shell.v !== this.view) return;
    if (this.currentModelRef == null) return;
    if (acc.getPane(EXPORT_PANE_NAME)) return;  // idempotent

    // Capture the ref AT CONSTRUCTION TIME of the pane factory, but read it lazily
    // so the dialog always sees the latest model text even if the pane sits around.
    const ownerRef = () => this.currentModelRef!;
    acc.addPane(EXPORT_PANE_NAME, () => renderExportPane(ownerRef), /* expanded = */ false);
  }),
);
```

Note the lazy `ownerRef`: the pane is added to an accordion once, but the accordion can be re-rendered and the user can keep editing the IVP text. By looking up `this.currentModelRef` *inside* the button handlers, the dialog always operates on the freshest text — no stale snapshot captured at pane-creation time.

`renderExportPane` takes a getter, not a value:

```typescript
// src/export/export-dialog.ts
export function renderExportPane(
  getRef: () => DiffStudioModelRef,
): HTMLElement {
  const open = (format: 'markdown' | 'latex') => {
    const ref = getRef();
    showExportDialog(ref.ivpText, ref.modelName, format);
  };
  const mdBtn = ui.button('Export as Markdown…', () => open('markdown'));
  const texBtn = ui.button('Export as LaTeX…', () => open('latex'));
  ui.tooltip.bind(mdBtn, 'Open the export dialog pre-configured for Markdown output.');
  ui.tooltip.bind(texBtn, 'Open the export dialog pre-configured for LaTeX output.');
  return ui.divV([mdBtn, texBtn], {style: {gap: '6px'}});
}
```

### 6.3 Teardown & scoping guarantees

Three layers of containment ensure the pane never leaks outside the Diff Studio session:

1. **Active-view guard** — `grok.shell.v !== this.view` early-returns. Accordions built while the user is in a different tab do not get the pane.
2. **Model-loaded guard** — `this.currentModelRef == null` early-returns. Diff Studio with a cold editor (no model yet) shows no Export pane.
3. **View-scoped lifetime** — the subscription lives in `this.view.subs`. `ViewBase.detach()` in [js-api/src/views/view.ts](../../../js-api/src/views/view.ts) iterates `subs` and unsubscribes each on close. Even if Datagrok caches accordions across activations, no *new* accordion built after Diff Studio closes can get the pane.
4. **Idempotency** — `acc.getPane(EXPORT_PANE_NAME)` prevents double-adding if the same accordion is re-constructed while Diff Studio is still active.

Caveat: if the platform **caches** an accordion that was built while Diff Studio was active and re-shows it later (same `acc.context`, same instance), the Export pane from the cached DOM may persist into the next activation of that object. This is worth verifying during smoke testing — if observed, we will add a `grok.events.onCurrentViewChanged` handler that strips `EXPORT_PANE_NAME` from any accordion whose construction preceded the view switch. Not expected in v1.

### 6.4 What does NOT change

- `src/package.ts` — untouched.
- The existing `@grok.decorators.fileHandler({ext: 'ipv'})` double-click behavior — untouched.
- The ribbon Download icon — untouched.
- No new package-level functions, no new `//tags: init`, no new `ObjectHandler`.

## Stage 7. Edge cases

**diff-grok throws.** Malformed `.ivp`: preview shows `// Error: <message>`, Copy and Download are disabled via `button.disabled = true` with a bound tooltip "Nothing to export — fix the errors in the model first."

**All Include checkboxes off.** If the result is structurally empty (no equations rendered), treat the same as a throw: disable actions with a tooltip.

**Model has no `#name`.** `modelName` falls back to the file base name or `'diff-export'` when assembling the default file name.

**Model has more than three equations and `compact` is on by default.** Compact mode still works but is less appropriate. The tooltip on **Compact layout** is explicit about the ≤3-equations rule; no auto-downgrade in v1.

**Large files.** `convertIvpToLatex` on the provided 13-example corpus runs in a few milliseconds. If the pollution model (20 species / 25 reactions) ever shows visible lag on regeneration, add a 50 ms debounce. Not expected in v1.

**Standalone toggle while format is Markdown.** `standaloneInput.root` is hidden via `display: none` whenever `format === 'markdown'`. Its value is preserved but ignored (`regenerate()` gates wrapping on `format === 'latex'`). Switching back to LaTeX restores the input with its previous state intact.

**Standalone + empty body.** If `convertIvpToLatex` returns an empty string (all Include boxes off), `wrapTex('')` still produces a syntactically valid — but empty — `.tex` document. Actions are disabled anyway by the "structurally empty" rule above, so the user cannot download/copy this artefact.

## Stage 8. Testing

- **Library side:** unchanged. All existing `diff-grok` Jest tests continue to cover conversion correctness.
- **Dialog side:** manual smoke test against every `.ivp` in `diff-grok/examples/`:
  - open dialog in both formats,
  - toggle every option, observe preview updates,
  - Copy, paste into a fresh Overleaf project (LaTeX) or Obsidian note (Markdown), confirm it renders,
  - Download, open the file, confirm LaTeX compiles with `pdflatex` out of the box and Markdown renders on GitHub.
- **Context-panel side:**
  - Open Diff Studio *without* a loaded model (cold editor) — confirm the **Diff Studio Export** pane does NOT appear anywhere in the Context panel.
  - Load a model (template, library example, user `.ivp`) — confirm the pane appears on the Context panel regardless of what `grok.shell.o` is (click the model itself, click a graph inside Diff Studio, click an output column — pane should be present in every accordion built while Diff Studio is active).
  - Switch to another app/view — confirm no Export pane shows up on newly constructed accordions there.
  - Switch back to Diff Studio — confirm the pane comes back on the next accordion.
  - Edit the IVP text in CodeMirror, then open the dialog from the pane — confirm the dialog's preview reflects the latest text (validates the lazy `getRef()` pattern).
  - Close the Diff Studio view, then in another app click on things — confirm the pane is gone (validates `view.subs` teardown).

## Order of work

1. Create `src/export/export-dialog.ts` skeleton — dialog, split, empty preview (Stage 2).
2. Add inputs and tooltips; wire `collectOptions()` (Stage 3).
3. Hook up live regeneration via `onChanged` (Stage 4).
4. Implement Copy and Download with local `wrapTex()` for LaTeX (Stage 5).
5. Create `DiffStudioModelRef` + `renderExportPane(getRef)`; maintain `this.currentModelRef` in the `DiffStudio` app class at every model-load hook. Wire `grok.events.onAccordionConstructed` from inside the app class with `grok.shell.v === this.view` + `this.currentModelRef != null` guard, and push the subscription onto `this.view.subs`. (Stage 6).
6. Handle error and empty states (Stage 7).
7. Smoke-test against the 13 examples (Stage 8).

Total new code: roughly 300–400 lines of TypeScript, split between `src/export/export-dialog.ts` (dialog + `wrapTex` + `renderExportPane`) and a small `DiffStudioModelRef` value class. The `DiffStudio` app class gains one subscription block (~6 lines) and one `grok.shell.o = …` assignment at the model-load hook. `src/package.ts` is not touched.

## Open points (answer during implementation)

- Exact Datagrok API for setting a tooltip on a button inside `ui.dialog` — use `ui.tooltip.bind()` if the plain `title` attribute is stripped by platform styling.
- Whether `ui.dialog` supports renaming the default OK button, or whether all three buttons (Cancel, Download, Copy) should be added via `addButton` for full control over label and order.
- Visual parity check: the preview `<pre>` must use the same monospace font token as other code blocks in Datagrok — confirm the CSS variable name during implementation.
