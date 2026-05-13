/**
 * Context-panel widget that renders the sense and antisense strands of an
 * OligoNucleotide cell as separate full molecular structures via Bio's
 * `Bio:toAtomicLevelPanel` function.
 *
 * The duplex HELM (e.g. `RNA1{...}|RNA2{...}$$$$`) is split into two
 * standalone HELM strings (`RNA1{...}$$$$` each), placed into a tiny
 * temporary DataFrame with proper Macromolecule/HELM tags, and each row's
 * cell is wrapped in a `DG.SemanticValue` and forwarded to `toAtomicLevelPanel`.
 *
 * Each strand is rendered as a collapsible accordion pane with lazy content:
 * the molfile assembly only runs when the pane is first expanded, so opening
 * a row's context panel never blocks on rendering antisense if the user only
 * cares about sense (and vice versa).
 */

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {canonicalizeHelm} from './helm-parser';

const HELM_SECTION_RE = /^(RNA|DNA|PEPTIDE|CHEM|BLOB)\d+\{/;

export function buildOligoStructuresPanel(value: DG.SemanticValue): DG.Widget {
  const helm: string = value.value ?? '';
  const polymerSection = helm.split('$')[0];
  const chains = polymerSection.split('|').map((c) => c.trim()).filter((c) => c.length > 0);

  if (chains.length === 0) {
    return DG.Widget.fromRoot(ui.divText('No HELM chains found', {
      style: {fontSize: '12px', color: 'var(--grey-4, #888)'},
    }));
  }

  // Renumber each chain to "RNA1{...}" and canonicalize aliased symbols
  // (e.g. `mR` → `m`, `fR` → `fl2r`, `LR` → `lna`, `sP` → `sp`) so the
  // central Bio monomer library — which only knows the canonical HELMCore
  // forms — can resolve every monomer when assembling the molfile.
  const standaloneHelms = chains.map((c) =>
    canonicalizeHelm(c.replace(HELM_SECTION_RE, '$11{') + '$$$$'),
  );
  const labels = chains.length === 1 ? ['Strand'] : ['Sense', 'Antisense'];

  // Build a temp DataFrame whose helm column carries the right tags so
  // Bio:toAtomicLevelPanel sees a real Macromolecule cell.
  const helmCol = DG.Column.fromStrings('helm', standaloneHelms);
  helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
  helmCol.meta.units = 'helm';
  helmCol.setTag('aligned', 'SEQ');
  helmCol.setTag('alphabet', 'RNA');
  helmCol.setTag('cell.renderer', 'helm');
  const df = DG.DataFrame.fromColumns([helmCol]);

  const acc = ui.accordion('oligo chemical structures');
  for (let i = 0; i < df.rowCount; i++) {
    const label = labels[i];
    const cell = df.cell(i, 'helm');
    // Lazy content: the molfile assembly fires only when the pane expands.
    acc.addPane(label, () => buildStrandPaneContent(cell, label));
  }
  return DG.Widget.fromRoot(acc.root);
}

/** Builds the host element for one strand's molecular structure. The actual
 * Bio:toAtomicLevelPanel call is async — the host is returned immediately so
 * the accordion expands without blocking, then content swaps in when ready. */
function buildStrandPaneContent(cell: DG.Cell, label: string): HTMLElement {
  const host = ui.div([], {style: {minHeight: '60px'}});
  host.appendChild(ui.divText('Building structure…', {
    style: {fontSize: '11px', color: 'var(--grey-4, #888)', padding: '4px'},
  }));

  (async () => {
    try {
      const sv = DG.SemanticValue.fromTableCell(cell);
      const widget = await grok.functions.call('Bio:toAtomicLevelPanel', {sequence: sv}) as DG.Widget;
      host.innerHTML = '';
      host.appendChild(widget?.root ?? notAvailable());
    } catch (e) {
      host.innerHTML = '';
      host.appendChild(ui.divText(
        `${label} render failed: ${(e as Error)?.message ?? String(e)}`,
        {style: {fontSize: '11px', color: 'var(--grey-4, #888)'}},
      ));
    }
  })();

  return host;
}

function notAvailable(): HTMLElement {
  return ui.divText('Structure not available', {
    style: {fontSize: '11px', color: 'var(--grey-4, #888)'},
  });
}
