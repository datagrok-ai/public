import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {splitAlignedPeptides} from '../utils/split-aligned';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Peptides} from '../peptides';
import {describe} from '../describe';
import {analyzePeptidesWidget} from '../widgets/analyze-peptides';
import {manualAlignmentWidget} from '../widgets/manual-alignment';
import {peptideMoleculeWidget} from '../widgets/peptide-molecule';
import * as P from '../package';

// let _package = new DG.Package();

category('peptides', async () => {
  let peptidesDf: DG.DataFrame;
  let options: StringDictionary;
  let peptidesGrid: DG.Grid;
  let asCol: DG.Column;
  let pepView: DG.TableView;

  before(async () => {
    // peptidesDf = DG.DataFrame.fromCsv(await P._package.files.readAsText('aligned.csv'));
    const csv = `ID,AlignedSequence,IC50
    1,NH2--A-Q-T-T-Y-K-N-Y-R-R-N-L-L--COOH,4.6411368455908086e-4
    2,NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH,0.003327324930165897
    3,NH2--A-N-T-T-Y-K-C-Y-R-R-N-L-L--COOH,3.0748148478921324e-4
    4,NH2--A-N-T-T-Y-K-F-Y-R-R-N-L-L--COOH,0.0015532837750281958
    5,NH2--A-V-T-T-Y-K-N-Y-R-R-N-L-L--COOH,6.549885174778741e-4
    6,NH2--A-N-T-T-Y-K-N-Y-R-R-N-L-L--COOH,0.00213298315038382
    7,NH2--A-N-T-T-Y-K-N-Y-R-F-N-L-L--COOH,0.002171297321903189
    8,NH2--A-N-T-T-Y-K-N-Y-R-R-N-H-L--COOH,0.002060711496394637
    9,NH2-M-A-N-T-T-Y-K-N-Y-R-R-N-L-L--COOH,0.0016058870359321516
    10,NH2--A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH,0.00212911675087592
    11,NH2--A-N-T-T-Y-K-N-Y-R-R-N-L-L--COOH,0.002736311013579287
    12,NH2--A-N-T-T-Y-K-N-Y-R-R-N-L-L--COOH,5.673074652436946e-5
    13,NH2-C-A-N-T-T-Y-K-N-Y-R-R-N-L-L--COOH,0.0032881139376902814
    14,NH2--A-N-T-T-Y-K-N-Y-R-H-N-L-L--COOH,0.0012828163841736553
    15,NH2-Y-A-N-T-T-Y-K-N-Y-R-D-N-L-L--COOH,7.186983807098166e-4
    16,NH2-M-A-N-T-T-Y-K-N-Y-R-N-N-L-L--COOH,0.00659708587488309
    17,NH2-P-A-N-T-T-Y-K-N-Y-R-G-N-L-L--COOH,3.7620528849324097e-4
    18,NH2-Y-A-N-T--Y-K-N-Y-R-S-N-L-L--COOH,6.812868474160967e-4
    19,NH2--A-N-T-T-Y-K-N-Y-R-S-N-L-L--COOH,0.0010148578953195436`;
    peptidesDf = DG.DataFrame.fromCsv(csv);
    options = {
      'activityColumnName': 'IC50',
      'scaling': '-lg',
    };
    asCol = peptidesDf.getCol('AlignedSequence');
    pepView = grok.shell.addTableView(peptidesDf);
    peptidesGrid = pepView.grid;
  });

  test('utils.split-sequence', async () => {
    splitAlignedPeptides(peptidesDf.getCol('AlignedSequence'));
  });

  test('describe', async () => {
    await describe(
      peptidesDf, options['activityColumnName'], options['scaling'], peptidesGrid, true,
      DG.BitSet.create(peptidesDf.rowCount, (i) => i % 2 === 0), true);
  });

  test('Peptides-class', async () => {
    const peptides = new Peptides();
    peptides.init(peptidesGrid, pepView, peptidesDf, options, asCol);
  });

  test('panel.peptides', async () => {
    await P.peptidesPanel(asCol);
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

  test('widgets.molfile', async () => {
    await P.peptideMolfile2('');
  });

  test('renderers.aligned-sequence-cell', async () => {
    P.alignedSequenceCellRenderer();
  });

  test('renderers.amino-acids-cell', async () => {
    P.aminoAcidsCellRenderer();
  });

  test('viewers.logo-viewer', async () => {
    P.logov();
  });

  test('viewers.sar-viewer', async () => {
    P.sar();
  });

  test('viewers.sar-vertical-viewer', async () => {
    P.sarVertical();
  });

  test('viewers.stacked-barchart-viewer', async () => {
    P.stackedBarChart();
  });

  test('viewers.subst-viewer', async () => {
    P.subst();
  });

  after(async () => {
    pepView.close();
    grok.shell.closeTable(peptidesDf);
  });
});
