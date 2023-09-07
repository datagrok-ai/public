import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../package';
import {startAnalysis} from '../widgets/peptides';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {PeptidesModel} from '../model';

export async function macromoleculeSarFastaDemoUI(): Promise<void> {
  const demo = new DemoScript('Macromolecule SAR analysis', '');
  const alignedSequenceCol = 'AlignedSequence';
  const simpleActivityColName = 'IC50';
  const simpleTable = DG.DataFrame.fromCsv(await _package.files.readAsText('aligned.csv'));
  const simpleActivityCol = simpleTable.getCol(simpleActivityColName);
  const simpleAlignedSeqCol = simpleTable.getCol(alignedSequenceCol);
  simpleAlignedSeqCol.semType = DG.SEMTYPE.MACROMOLECULE;
  simpleAlignedSeqCol.setTag(C.TAGS.ALPHABET, ALPHABET.PT);
  simpleAlignedSeqCol.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
  simpleAlignedSeqCol.setTag(bioTAGS.aligned, ALIGNMENT.SEQ_MSA);
  const simpleScaledCol = scaleActivity(simpleActivityCol, C.SCALING_METHODS.MINUS_LG);

  demo.step('Load data', async () => {grok.shell.addTableView(simpleTable);},
    {description: 'Load the dataset containing macromolecule sequences. Notice how Datagrok detects macromolecules ' +
      'and applies renderer for better visual understanding of the data', delay: 2000});

  let alignedCol: DG.Column<string>;
  demo.step('Align sequences', async () => {
    const clustersCol = DG.Column.string('Cluster', simpleTable.rowCount).init('0');
    alignedCol = await grok.functions.call('Bio:alignSequences', {sequenceCol: simpleAlignedSeqCol, clustersCol});
  }, {description: 'Align sequences with our Multiple Sequence Alignment tool. This tool can be found in the top ' +
    'menu Bio -> Alignment -> MSA... New msa column will be added to the table', delay: 2000});

  let model: PeptidesModel;
  demo.step('Run SAR analysis', async () => {
    model = await startAnalysis(simpleActivityCol, alignedCol, null, simpleTable, simpleScaledCol,
      C.SCALING_METHODS.MINUS_LG) as PeptidesModel;
    model.analysisView.addViewer('WebLogo');
  }, {description: 'Run SAR analysis on aligned sequences from top menu: Bio -> SAR -> Peptides...', delay: 2000});

  demo.step('Browse Mutation Cliffs', async () => {model.modifyMutationCliffsSelection({positionOrClusterType: 'D', monomerOrCluster: '13'});},
    {description: 'Browse Mutation Cliffs by selecting corresponding cells in Mutation Cliffs viewer in the bottom',
      delay: 2000});
  await demo.start();
}
