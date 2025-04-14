import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_package, getSeqHelper, toAtomicLevel} from '../package';


export async function toAtomicLevelWidget(sequence: DG.SemanticValue): Promise<DG.Widget> {
  const errorDiv = ui.divText('');
  const errorWidget = DG.Widget.fromRoot(errorDiv);
  try {
    if (!sequence || !sequence.value) {
      errorDiv.innerText = 'No sequence provided';
      return errorWidget;
    }
    if (!sequence.cell || !sequence.cell.dart || !sequence.cell.dataFrame || !sequence.cell.column) {
      errorDiv.innerText = 'Atomic level conversion requeires a sequence column';
      return errorWidget;
    }
    const supportedUnits: string[] = [NOTATION.FASTA, NOTATION.SEPARATOR, NOTATION.HELM];
    //todo: add support for custom notations
    if (!supportedUnits.includes(sequence.cell.column.meta.units?.toLowerCase() ?? '')) {
      errorDiv.innerText = 'Unsupported sequence notation. please use Bio | Polytool | Convert';
      return errorWidget;
    }
    const seqHelper = await getSeqHelper();
    const seqSh = seqHelper.getSeqHandler(sequence.cell.column);
    if (!seqSh) {
      errorDiv.innerText = 'No sequence handler found';
      return errorWidget;
    }
    if ((seqSh.getSplitted(sequence.cell.rowIndex, 50)?.length ?? 100) > 40) {
      errorDiv.innerText = 'Maximum number of monomers is 40';
      return errorWidget;
    }
    const singleValCol = DG.Column.fromStrings('singleVal', [sequence.value]);
    const sDf = DG.DataFrame.fromColumns([singleValCol]);
    // copy over all the tags
    Object.entries(sequence.cell.column.tags).forEach(([key, value]) => {
      singleValCol.setTag(key, value as string);
    });
    await toAtomicLevel(sDf, singleValCol, sequence.cell.column.meta.units === NOTATION.HELM, false);
    if (sDf.columns.length < 2) {
      errorDiv.innerText = 'No structure generated';
      return errorWidget;
    }
    const molCol = sDf.columns.byIndex(1);
    const molfile = molCol.get(0);
    if (!molfile) {
      errorDiv.innerText = 'No structure generated';
      return errorWidget;
    }
    molCol.semType = DG.SEMTYPE.MOLECULE;
    const molSemanticValue = DG.SemanticValue.fromTableCell(sDf.cell(0, molCol.name));
    const panel = ui.panels.infoPanel(molSemanticValue);
    let molPanel: DG.Widget | null = null;
    if (panel)
      molPanel = DG.Widget.fromRoot(panel.root);


    const root = grok.chem.drawMolecule(molfile, 300, 300, false);
    root.style.cursor = 'pointer';
    ui.tooltip.bind(root, 'Click to expand');
    root.onclick = () => {
      const width = window.innerWidth - 200;
      const height = window.innerHeight - 200;
      const bigMol = grok.chem.drawMolecule(molfile, width, height, false);
      ui.dialog({title: 'Molecule'}).add(bigMol).showModal(true);
    };
    if (molPanel)
      molPanel.root.prepend(root);
    return molPanel ?? DG.Widget.fromRoot(root);
  } catch (e) {
    _package.logger.error(e);
  }

  errorDiv.innerText = 'No Structure generated';
  return errorWidget;
}
