/**
 * Shared structured-tooltip builder (`.nca-richtip`).
 *
 * A single tooltip look across the whole app — form inputs (mapping pane) AND
 * every grid header (results, aggregated, study, TK, bioavailability,
 * dose-proportionality). Long, multi-clause explanations read poorly as a flat
 * `tooltipText` / bare `ui.divText` run-together string; this builder renders
 * them with visual hierarchy: bold header → summary lede → optional formula /
 * detail / convention note / default, plus optional trailing `extra` rows (a
 * colour legend, typically). Styling lives in `.nca-richtip*` (css/app-styles.css).
 *
 * Extracted from `nca-mapping-pane.ts` (study-grid-header-tooltips) so the grids —
 * which call `ui.tooltip.show(richTipEl(...))` — and the inputs — which bind
 * `richTip(...)` via {@link bindInputTip} — share ONE renderer rather than each
 * hand-rolling a `ui.divV([ui.h2(title), ui.divText(desc)])` with no styling.
 */

import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';

/** Spec for a structured ("rich") tooltip. Render order: header → summary →
 *  formula → detail → note → default → extra; empty fields skipped. */
export interface RichTipSpec {
  /** Name / label — rendered as the bold header. */
  readonly title: string;
  /** One-line "what it is" — the lede. */
  readonly summary: string;
  /** Optional formula / expression, rendered monospace on its own line. */
  readonly formula?: string;
  /** Optional "why / when to use" paragraph (muted). */
  readonly detail?: string;
  /** Optional standard / convention reference (muted italic, e.g. CDISC / ICH),
   *  or an in-grid formatting rule (e.g. "Bold when n < 4"). */
  readonly note?: string;
  /** Optional default-value note, rendered as `Default: …` (muted, smallest). */
  readonly def?: string;
  /** Optional trailing rows appended verbatim after the text block — a colour
   *  legend (`legendSpacer()` + `legendRow()`s) for the conditional-format grids. */
  readonly extra?: ReadonlyArray<HTMLElement>;
}

/** Build the tooltip element. Grids call this directly for `ui.tooltip.show`. */
export function richTipEl(spec: RichTipSpec): HTMLElement {
  const rows: HTMLElement[] = [
    ui.h3(spec.title),
    ui.divText(spec.summary, 'nca-richtip-summary'),
  ];
  if (spec.formula !== undefined) rows.push(ui.divText(spec.formula, 'nca-richtip-formula'));
  if (spec.detail !== undefined) rows.push(ui.divText(spec.detail, 'nca-richtip-detail'));
  if (spec.note !== undefined) rows.push(ui.divText(spec.note, 'nca-richtip-note'));
  if (spec.def !== undefined) rows.push(ui.divText(`Default: ${spec.def}`, 'nca-richtip-default'));
  if (spec.extra !== undefined) rows.push(...spec.extra);
  return ui.divV(rows, 'nca-richtip');
}

/** Lazy factory for `ui.tooltip.bind` (inputs, section headers), which accepts
 *  only `string | () => Element`. */
export function richTip(spec: RichTipSpec): () => HTMLElement {
  return (): HTMLElement => richTipEl(spec);
}

/**
 * Bind a rich tooltip to a form input. DG inputs own a native tooltip channel
 * (`tooltipText` → `InputBase.setTooltip`, string-only) whose hover handler masks
 * a DOM `ui.tooltip.bind` on `input.root` for non-toggle inputs — which is why
 * tooltips silently vanish there. The platform instead honours DOM tooltips bound
 * to the input's own child elements (the pattern DG itself uses, e.g. compute-utils
 * `fitting-view`). Bind BOTH the caption label and the editor so the whole input
 * row is a hover target, matching the old native `tooltipText`.
 */
export function bindInputTip(input: DG.InputBase, tip: () => HTMLElement): void {
  ui.tooltip.bind(input.captionLabel, tip);
  ui.tooltip.bind(input.input, tip);
}
