// Context-panel handler for `PL Diagram` cells (per-row LigNetwork
// thumbnails written by `runPlBatch`).
//
// Registered once via `DG.ObjectHandler.register(new PlDiagramObjectHandler())`
// in BSV's `init()`. The platform routes cell clicks through every
// applicable handler's `renderProperties` and merges results into the
// right sidebar. Claims ownership of cells whose column carries the
// `%prolifSource` tag (value = name of the source ligand/PDB column,
// set by `runPlBatch`). Cell rendering itself is delegated to PowerGrid's
// built-in `rawPng` renderer via `semType: 'rawPng'` on the column.

import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getPlHtmlForRow, interactionsColForDiagram, renderInteractionBreakdown} from './prolif-panel';

/** Tag set on every `PL Diagram` column by the batch handler. Its value is
 *  the name of the source ligand/PDB column the diagrams were computed
 *  from. Presence of this tag is what `PlDiagramObjectHandler.isApplicable`
 *  uses to claim ownership of the cell. */
export const PROLIF_SOURCE_TAG = '%prolifSource';

export class PlDiagramObjectHandler extends DG.ObjectHandler {
  get type(): string { return 'pl-diagram'; }

  isApplicable(x: any): boolean {
    if (!(x instanceof DG.SemanticValue)) return false;
    const cell = x.cell;
    return cell != null && cell.dart != null && cell.column != null &&
      cell.column.tags[PROLIF_SOURCE_TAG] != null;
  }

  renderProperties(x: DG.SemanticValue, _context: unknown = null): HTMLElement {
    // `SemanticValue.cell` is the typical path; `grok.shell.o` covers the
    // case where the panel re-fires without a fresh cell reference (e.g.
    // user clicks elsewhere then back, and the SemanticValue is stale).
    const dgCell = x.cell ?? (grok.shell.o as DG.Cell | null);
    const df = dgCell?.dataFrame ?? grok.shell.t ?? null;
    const rowIdx = dgCell?.rowIndex ?? df?.currentRowIdx ?? -1;
    if (df == null || rowIdx < 0)
      return ui.divText('No row context.');
    const html = getPlHtmlForRow(df, rowIdx);
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
    const interactionsCol = dgCell?.column?.name != null
      ? interactionsColForDiagram(df, dgCell.column.name) : null;
    const interactionsStr = (interactionsCol?.get(rowIdx) as string | null) ?? '';
    return ui.div(
      [iframe, renderInteractionBreakdown(interactionsStr)], 'd4-empty-parent');
  }
}
