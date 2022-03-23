import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {PeptidesController} from '../peptides';
import {analyzePeptidesWidget} from '../widgets/analyze-peptides';
import {manualAlignmentWidget} from '../widgets/manual-alignment';
import {peptideMoleculeWidget} from '../widgets/peptide-molecule';
import {_packageTest} from '../package-test';

category('peptides', async () => {
  let peptidesDf: DG.DataFrame;
  let options: StringDictionary;
  let peptidesGrid: DG.Grid;
  let asCol: DG.Column;
  let pepView: DG.TableView;

  before(async () => {
    peptidesDf = DG.DataFrame.fromCsv(await _packageTest.files.readAsText('aligned.csv'));
    options = {
      activityColumnName: 'IC50',
      scaling: '-lg',
    };
    asCol = peptidesDf.getCol('AlignedSequence');
    pepView = grok.shell.addTableView(peptidesDf);
    peptidesGrid = pepView.grid;
  });

  test('utils.split-sequence', async () => {
    PeptidesController.splitAlignedPeptides(peptidesDf.getCol('AlignedSequence'));
  });

  // test('describe', async () => {
  //   await describe(
  //     peptidesDf, options['activityColumnName'], options['scaling'], peptidesGrid, true,
  //     DG.BitSet.create(peptidesDf.rowCount, (i) => i % 2 === 0), true);
  // });

  test('Peptides-controller', async () => {
    const peptides = await PeptidesController.getInstance(peptidesDf);
    peptides.init(peptidesGrid, pepView, options, asCol, peptidesDf.columns.names());
  });

  test('widgets.analyze-peptides', async () => {
    await analyzePeptidesWidget(asCol, pepView, peptidesGrid, peptidesDf);
  });

  test('widgets.manual-alignment', async () => {
    manualAlignmentWidget(asCol, peptidesDf);
  });

  test('widgets.peptide-molecule', async () => {
    await peptideMoleculeWidget('NH2--A-N-T-T-Y-K-N-Y-R-S-N-L-L--COOH');
  });

  after(async () => {
    pepView.close();
    grok.shell.closeTable(peptidesDf);
  });
});
