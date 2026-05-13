// Context-panel handler for `PL Diagram` cells (per-row LigNetwork
// thumbnails written by `runPlBatch`).
//
// Registered once via `DG.ObjectHandler.register(new PlDiagramObjectHandler())`
// in BSV's `init()`. The platform routes cell clicks through every
// applicable handler's `renderProperties` and merges results into the
// right sidebar. Claims ownership of cells whose column carries the
// `.%prolif-source` tag (value = name of the source ligand/PDB column,
// set by `runPlBatch`). Cell rendering itself is delegated to PowerGrid's
// built-in `rawPng` renderer via `semType: 'rawPng'` on the column.
//
// Tag naming: `.%`-prefixed kebab-case matches the Datagrok convention
// for column-linking tags (see Chem's `.%chem-scaffold-align`,
// `.%chem-space-embedding-col`, Curves' `.%curve-format`).

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  getPlHtmlForRow, interactionsColForDiagram,
  PL_DIAGRAM_SEM_TYPE, PROLIF_SOURCE_TAG, renderInteractionBreakdown,
} from './prolif';

// Re-export so existing call sites and tests keep working.
export {PROLIF_SOURCE_TAG};

export class PlDiagramObjectHandler extends DG.ObjectHandler {
  get type(): string { return 'pl-diagram'; }

  isApplicable(x: any): boolean {
    if (!(x instanceof DG.SemanticValue)) return false;
    const cell = x.cell;
    return cell != null && cell.dart != null && cell.column != null && x.semType === PL_DIAGRAM_SEM_TYPE && 
      cell.column.tags[PROLIF_SOURCE_TAG] != null && !!cell.dataFrame.col(cell.column.tags[PROLIF_SOURCE_TAG]);
  }

  renderProperties(x: DG.SemanticValue, _context: unknown = null): HTMLElement {
    const dgCell = x.cell;
    const diagramCol = dgCell?.column ?? null;
    const df = dgCell?.dataFrame ?? null;
    const rowIdx = dgCell?.rowIndex ?? -1;
    if (df == null || diagramCol == null || rowIdx < 0)
      return ui.divText('No row context.');
    // HTML lives on the diagram column's `temp` slot (populated by
    // `runPlBatch`). Survives until the column is removed or the
    // DataFrame is closed; doesn't survive project save/reopen.
    const html = getPlHtmlForRow(diagramCol, rowIdx);
    if (!html) {
      return ui.divText(
        `LigNetwork HTML not cached for ${df.name} row ${rowIdx} — ` +
        `re-run "Compute for whole dataset" to refresh.`,
      );
    }
    const iframe = ui.element('iframe') as HTMLIFrameElement;
    iframe.srcdoc = html;
    iframe.classList.add('bsv-pl-panel-iframe');
    iframe.onload = () => iframe.classList.add('bsv-pl-panel-iframe-loaded');
    const interactionsCol = interactionsColForDiagram(diagramCol);
    const interactionsStr = (interactionsCol?.get(rowIdx) as string | null) ?? '';
    const originalAcc = ui.panels.infoPanel(x);
    originalAcc.addPane('Interactions', () => {
      return ui.div(
        [iframe, renderInteractionBreakdown(interactionsStr)], 'd4-empty-parent');
    }, true);
    return originalAcc.root;
  }
}
