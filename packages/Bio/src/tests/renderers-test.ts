import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {category, expect, test, awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

import {importFasta} from '../package';
import {convertDo} from '../utils/convert';
import {generateLongSequence, generateManySequences, performanceTest} from './utils/sequences-generators';
import {multipleSequenceAlignmentUI} from '../utils/multiple-sequence-alignment-ui';
import {awaitGrid} from './utils';
import * as C from '../utils/constants';

category('renderers', () => {
  test('long sequence performance ', async () => {
    await performanceTest(generateLongSequence, 'Long sequences');
  });

  test('many sequence performance', async () => {
    await performanceTest(generateManySequences, 'Many sequences');
  });
  test('many sequence performance', async () => {
    await performanceTest(generateManySequences, 'Many sequences');
  });

  test('rendererMacromoleculeFasta', async () => {
    await _rendererMacromoleculeFasta();
  });

  test('rendererMacromoleculeSeparator', async () => {
    await _rendererMacromoleculeSeparator();
  });

  test('rendererMacromoleculeDifference', async () => {
    await _rendererMacromoleculeDifference();
  });

  test('afterMsa', async () => {
    await _testAfterMsa();
  });

  test('afterConvert', async () => {
    await _testAfterConvert();
  });

  test('selectRendererBySemType', async () => {
    await _selectRendererBySemType();
  });

  async function _rendererMacromoleculeFasta() {
    const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/FASTA.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);

    const seqCol = df.getCol('Sequence');
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    const resCellRenderer = seqCol.getTag(DG.TAGS.CELL_RENDERER);
    expect(resCellRenderer, 'sequence');
  }

  async function _rendererMacromoleculeSeparator() {
    const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/SEPARATOR_PT.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);

    const seqCol = df.getCol('sequence');
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: seqCol});
    if (semType)
      seqCol.semType = semType;

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    const resCellRenderer = seqCol.getTag(DG.TAGS.CELL_RENDERER);
    expect(resCellRenderer, 'sequence');
  }

  async function _rendererMacromoleculeDifference() {
    const seqDiffCol: DG.Column = DG.Column.fromStrings('SequencesDiff',
      ['meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca#D-Tyr_Et/Tyr_ab-dehydroMe/meN/E/N/dV']);
    seqDiffCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    seqDiffCol.setTag(bioTAGS.separator, '/');
    seqDiffCol.setTag(bioTAGS.aligned, 'SEQ');
    seqDiffCol.setTag(bioTAGS.alphabet, 'UN');
    seqDiffCol.setTag(bioTAGS.alphabetIsMultichar, 'true');
    seqDiffCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
    const df = DG.DataFrame.fromColumns([seqDiffCol]);

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    const resCellRenderer = seqDiffCol.getTag(DG.TAGS.CELL_RENDERER);
    expect(resCellRenderer, C.SEM_TYPES.MACROMOLECULE_DIFFERENCE);
  }

  async function _testAfterMsa() {
    const fastaTxt: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/FASTA.fasta');
    const df: DG.DataFrame = importFasta(fastaTxt)[0];

    const srcSeqCol: DG.Column = df.getCol('sequence');
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: srcSeqCol});
    if (semType)
      srcSeqCol.semType = semType;

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    console.log('Bio: tests/renderers/afterMsa, table view');
    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    console.log('Bio: tests/renderers/afterMsa, src before test ' +
      `semType="${srcSeqCol!.semType}", units="${srcSeqCol!.getTag(DG.TAGS.UNITS)}", ` +
      `cell.renderer="${srcSeqCol!.getTag(DG.TAGS.CELL_RENDERER)}"`);
    expect(srcSeqCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(srcSeqCol.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(srcSeqCol.getTag(bioTAGS.aligned), ALIGNMENT.SEQ);
    expect(srcSeqCol.getTag(bioTAGS.alphabet), ALPHABET.PT);
    expect(srcSeqCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    const msaSeqCol = await multipleSequenceAlignmentUI({col: srcSeqCol});
    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    expect(msaSeqCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(msaSeqCol.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(msaSeqCol.getTag(bioTAGS.aligned), ALIGNMENT.SEQ_MSA);
    expect(msaSeqCol.getTag(bioTAGS.alphabet), ALPHABET.PT);
    expect(msaSeqCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    // check newColumn with UnitsHandler constructor
    const _uh: UnitsHandler = UnitsHandler.getOrCreate(msaSeqCol);
  }

  async function _testAfterConvert() {
    const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/FASTA_PT.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);

    const srcCol: DG.Column = df.col('sequence')!;
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: srcCol});
    if (semType)
      srcCol.semType = semType;

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    const tgtCol: DG.Column = await convertDo(srcCol, NOTATION.SEPARATOR, '/');
    await awaitGrid(tv.grid);
    expect(tv.grid.dataFrame.id, df.id);

    const resCellRenderer = tgtCol.getTag(DG.TAGS.CELL_RENDERER);
    expect(resCellRenderer, 'sequence');

    // check tgtCol with UnitsHandler constructor
    const _uh: UnitsHandler = UnitsHandler.getOrCreate(tgtCol);
  }

  async function _selectRendererBySemType() {
    /* There are renderers for semType Macromolecule and MacromoleculeDifference.
       Misbehavior was by selecting Macromolecule renderers for MacromoleculeDifference semType column
    /**/
    const seqDiffCol: DG.Column = DG.Column.fromStrings('SequencesDiff',
      ['meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca#D-Tyr_Et/Tyr_ab-dehydroMe/meN/E/N/dV']);
    seqDiffCol.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    seqDiffCol.setTag(bioTAGS.separator, '/');
    seqDiffCol.setTag(bioTAGS.aligned, 'SEQ');
    seqDiffCol.setTag(bioTAGS.alphabet, 'UN');
    seqDiffCol.setTag(bioTAGS.alphabetIsMultichar, 'true');
    seqDiffCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
    const df = DG.DataFrame.fromColumns([seqDiffCol]);
    const tv = grok.shell.addTableView(df);

    await delay(100);
    const renderer = seqDiffCol.getTag(DG.TAGS.CELL_RENDERER);
    if (renderer !== 'MacromoleculeDifference') { // this is value of MacromoleculeDifferenceCR.cellType
      throw new Error(`Units 'separator', separator '/' and semType 'MacromoleculeDifference' ` +
        `have been manually set on column but after df was added as table, ` +
        `view renderer has set to '${renderer}' instead of correct 'MacromoleculeDifference'.`);
    }
  }
});
