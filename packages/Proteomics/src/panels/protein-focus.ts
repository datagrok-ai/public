import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';

const PROTEIN_ID_HINTS =
  ['primary protein id', 'pg.proteingroups', 'protein id', 'uniprot', 'accession'];

/**
 * Selects a protein row so the semType-registered UniProt panel populates in the
 * context panel, opens the context panel, and (best-effort) expands the UniProt
 * pane so the result is actually visible rather than one click away.
 *
 * Shared by two callers: every import handler passes the default row 0 (a small
 * "it works" cue the moment a file opens), and the demo passes its top
 * differential-expression hit. No-op if the table has no recognizable protein-id
 * column, so it's safe to call after any import.
 */
export function focusProtein(df: DG.DataFrame, rowIdx: number = 0): void {
  if (df.rowCount === 0 || rowIdx < 0 || rowIdx >= df.rowCount) return;
  const idCol = findColumn(df, SEMTYPE.PROTEIN_ID, PROTEIN_ID_HINTS);
  if (!idCol) return;
  df.currentCell = df.cell(rowIdx, idCol.name);
  grok.shell.windows.showContextPanel = true;
  void expandUniProtPane();
}

/**
 * Best-effort: expands the UniProt widget in the semantic context panel so the
 * protein result is visible immediately.
 *
 * There is no supported API to expand a semantic sub-accordion pane from a
 * package function — `onAccordionConstructed` only surfaces the top-level
 * `Actions` accordion, while the UniProt widget is rendered inside the platform's
 * semantic panel — so this reaches into the context-panel DOM. It is EXPAND-ONLY
 * (it clicks a pane header only when it lacks the `expanded` class, so it never
 * collapses a pane the user already opened) and fully guarded: if the platform
 * DOM ever changes shape, the panel simply stays collapsed, exactly as if this
 * helper were absent. The panel content itself is populated either way by the
 * current-cell selection in `focusProtein`.
 */
async function expandUniProtPane(): Promise<void> {
  const wanted = ['Proteomics', 'UniProt']; // group pane, then the widget inside it
  for (let attempt = 0; attempt < 12; attempt++) {
    await new Promise((r) => setTimeout(r, 250));
    try {
      const acc = document.querySelector('[d4-title="semantic meta.Proteomics-ProteinId"]');
      if (!acc) continue;
      let allExpanded = true;
      for (const name of wanted) {
        const header = Array.from(acc.querySelectorAll('.d4-accordion-pane-header'))
          .find((h) => h.textContent?.trim() === name) as HTMLElement | undefined;
        if (!header) { allExpanded = false; continue; }
        if (!header.classList.contains('expanded')) {
          header.click(); // reveals the child pane on a later iteration
          allExpanded = false;
        }
      }
      if (allExpanded) return;
    } catch { /* best-effort presentation only — leave collapsed on any DOM change */ }
  }
}
