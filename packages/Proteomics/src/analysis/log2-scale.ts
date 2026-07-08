import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';
import {log2TransformColumns, copyAsLog2Columns, detectLog2Status} from '../parsers/shared-utils';

/**
 * Post-import correction for the log2 scale of intensity data.
 *
 * At import each parser decides whether the intensities are raw (needs a log2
 * transform) or already log2-transformed (copy as-is). That decision is a
 * magnitude heuristic (`detectLog2Status`) and has misfired on small-magnitude
 * linear intensities — values that happen to sit in the [0, 30] "looks like
 * log2" band get copied instead of transformed, producing nonsense ±thousands
 * fold changes downstream. This module lets the analyst override that call
 * without re-importing.
 *
 * The rebuild works off the ORIGINAL (non-`log2(...)`) intensity columns, which
 * every parser keeps pristine — normalization/imputation only ever mutate the
 * `log2(...)` copies in place. So re-deriving the `log2(...)` columns is always
 * valid and, by design, discards any normalization/imputation/DE that ran on the
 * old values; those steps must be re-run, which the dialog warns about.
 */

/** Original (pre-`log2(...)`) intensity column names — the pristine source the
 * rebuild derives from. Excludes the `log2(...)` copies themselves. */
export function getIntensityOriginals(df: DG.DataFrame): string[] {
  return df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && !c.name.startsWith('log2('))
    .map((c) => c.name);
}

/** Infers whether the current `log2(...)` columns were produced by a transform
 * (`Math.log2(orig)`) or a straight copy (`orig`, i.e. treated as already log2).
 * Scans until the first decisive non-null positive row — the two branches never
 * coincide for a positive value (`log2(x) ≠ x` for all x > 0). Returns 'unknown'
 * when no `log2(...)` column or no decisive row exists. */
export function detectCurrentLog2Applied(
  df: DG.DataFrame, originalNames: string[],
): 'transform' | 'copy' | 'unknown' {
  for (const name of originalNames) {
    const orig = df.col(name);
    const l2 = df.col(`log2(${name})`);
    if (!orig || !l2) continue;
    for (let i = 0; i < df.rowCount; i++) {
      if (orig.isNone(i) || l2.isNone(i)) continue;
      const o = Number(orig.get(i));
      const v = Number(l2.get(i));
      if (!(o > 0) || !Number.isFinite(v)) continue;
      if (Math.abs(v - Math.log2(o)) < 1e-4) return 'transform';
      if (Math.abs(v - o) < 1e-4) return 'copy';
    }
  }
  return 'unknown';
}

/** Rebuilds the `log2(...)` columns from the originals under the chosen scale:
 * `alreadyLog2` copies as-is, otherwise applies the log2 transform. Removes the
 * existing `log2(...)` columns first (so a re-run doesn't duplicate) and clears
 * stale downstream pipeline tags — the base data changed, so any prior
 * normalize/impute/DE is invalid and must be re-run. No-op (returns false) when
 * the requested scale already matches what's applied. */
export function applyLog2Scale(df: DG.DataFrame, alreadyLog2: boolean): boolean {
  const originals = getIntensityOriginals(df);
  if (originals.length === 0) return false;

  const current = detectCurrentLog2Applied(df, originals);
  const target = alreadyLog2 ? 'copy' : 'transform';
  if (current === target) return false;

  for (const name of originals) {
    const l2 = `log2(${name})`;
    if (df.col(l2)) df.columns.remove(l2);
  }
  if (alreadyLog2)
    copyAsLog2Columns(df, originals);
  else
    log2TransformColumns(df, originals);

  // Base intensities were re-derived — anything computed from the old log2
  // columns is stale. Drop the completion tags so viewers gate correctly and
  // the analyst re-runs from Normalize.
  for (const tag of ['proteomics.normalized', 'proteomics.imputed',
    'proteomics.de_complete', 'proteomics.de_method'])
    df.setTag(tag, '');
  return true;
}

/** Dialog to review and override the log2 scale of the imported intensities.
 * Seeds the checkbox from the scale actually applied at import (inferred), not a
 * re-detection, so it reflects reality; shows the magnitude heuristic's read as
 * an advisory hint. */
export function showLog2ScaleDialog(df: DG.DataFrame): void {
  const originals = getIntensityOriginals(df);
  if (originals.length === 0) {
    grok.shell.warning('No intensity columns found — nothing to rescale.');
    return;
  }

  const current = detectCurrentLog2Applied(df, originals);
  const detection = detectLog2Status(df, originals);

  const alreadyLog2 = ui.input.bool('Data is already log2-transformed', {
    value: current === 'copy',
    tooltipText: 'On: intensities are used as-is. Off: a log2 transform is applied. ' +
      'Flip this if fold changes look wildly off — the import heuristic can misread ' +
      'small-magnitude raw intensities as already-log2.',
  });

  const hint = ui.divText(`Auto-detection: ${detection.message}.` +
    (current !== 'unknown' ? ` Currently applied: ${current === 'copy' ?
      'used as-is (no transform)' : 'log2-transformed'}.` : ''));
  hint.style.cssText = 'font-style:italic;color:#888;margin:4px 0;font-size:12px;max-width:420px;';

  const pipelineRan = ['proteomics.normalized', 'proteomics.imputed', 'proteomics.de_complete']
    .some((t) => df.getTag(t) === 'true');
  const warn = ui.divText(
    'Normalize / Impute / Differential Expression have already run. Changing the ' +
    'scale resets those steps — re-run them from Analyze.');
  warn.style.cssText = 'color:#b26a00;margin:4px 0;font-size:12px;max-width:420px;';
  warn.style.display = pipelineRan ? '' : 'none';

  ui.dialog('Set Log2 Scale')
    .add(alreadyLog2)
    .add(hint)
    .add(warn)
    .onOK(() => {
      const changed = applyLog2Scale(df, alreadyLog2.value);
      if (!changed) {
        grok.shell.info('Scale unchanged.');
        return;
      }
      grok.shell.info(alreadyLog2.value ?
        'Intensities are now used as-is (no log2 transform). Re-run analysis steps.' :
        'Intensities are now log2-transformed. Re-run analysis steps.');
    })
    .show();
}
