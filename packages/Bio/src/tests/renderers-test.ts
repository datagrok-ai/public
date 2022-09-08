import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {importFasta, multipleSequenceAlignmentAny} from '../package';
import {convertDo} from '../utils/convert';
import {ALPHABET, NOTATION, UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {SEM_TYPES, TAGS} from '../utils/constants';

category('renderers', () => {
  let tvList: DG.TableView[];
  let dfList: DG.DataFrame[];

  before(async () => {
    await grok.functions.call('Bio:initBio');
    tvList = [];
    dfList = [];
  });

  after(async () => {
    dfList.forEach((df: DG.DataFrame) => { grok.shell.closeTable(df); });
    tvList.forEach((tv: DG.TableView) => tv.close());
  });

  test('afterMsa', async () => {
    await _testAfterMsa();
  });

  test('afterConvert', async () => {
    await _testAfterConvert();
  });

  test('setRenderer', async () => {
    await _setRendererManually();
  });

  async function _testAfterMsa() {
    const fastaTxt: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.fasta');
    const df: DG.DataFrame = importFasta(fastaTxt)[0];

    const srcSeqCol: DG.Column = df.getCol('sequence');
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: srcSeqCol});
    if (semType)
      srcSeqCol.semType = semType;

    const tv: DG.TableView = grok.shell.addTableView(df);
    // call to calculate 'cell.renderer' tag
    await grok.data.detectSemanticTypes(df);

    console.log('Bio: tests/renderers/afterMsa, table view');

    console.log('Bio: tests/renderers/afterMsa, src before test ' +
      `semType="${srcSeqCol!.semType}", units="${srcSeqCol!.getTag(DG.TAGS.UNITS)}", ` +
      `cell.renderer="${srcSeqCol!.getTag(DG.TAGS.CELL_RENDERER)}"`);
    expect(srcSeqCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(srcSeqCol.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(srcSeqCol.getTag(UnitsHandler.TAGS.aligned), 'SEQ');
    expect(srcSeqCol.getTag(UnitsHandler.TAGS.alphabet), ALPHABET.PT);
    expect(srcSeqCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    const msaSeqCol: DG.Column = (await multipleSequenceAlignmentAny(df, srcSeqCol!))!;
    tv.grid.invalidate();

    expect(msaSeqCol.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(msaSeqCol.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(msaSeqCol.getTag(UnitsHandler.TAGS.aligned), 'SEQ.MSA');
    expect(msaSeqCol.getTag(UnitsHandler.TAGS.alphabet), ALPHABET.PT);
    expect(msaSeqCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    // check newColumn with UnitsHandler constructor
    const uh: UnitsHandler = new UnitsHandler(msaSeqCol);

    dfList.push(df);
    tvList.push(tv);
  }

  async function _testAfterConvert() {
    const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA_PT.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const tv: DG.TableView = grok.shell.addTableView(df);

    const srcCol: DG.Column = df.col('sequence')!;
    // await grok.data.detectSemanticTypes(df);
    const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: srcCol});
    if (semType)
      srcCol.semType = semType;
    await grok.data.detectSemanticTypes(df);

    const tgtCol: DG.Column = await convertDo(srcCol, NOTATION.SEPARATOR, '/');
    expect(tgtCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    // check tgtCol with UnitsHandler constructor
    const uh: UnitsHandler = new UnitsHandler(tgtCol);

    tvList.push(tv);
    dfList.push(df);
  }

  async function _setRendererManually() {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings(
      'SequencesDiff', ['meI/hHis/Aca/N/T/dK/Thr_PO3H2/Aca#D-Tyr_Et/Tyr_ab-dehydroMe/meN/E/N/dV'])]);
    df.col('SequencesDiff')!.tags[DG.TAGS.UNITS] = 'separator';
    df.col('SequencesDiff')!.tags[TAGS.SEPARATOR] = '/';
    df.col('SequencesDiff')!.semType = SEM_TYPES.MACROMOLECULE_DIFFERENCE;
    const tw = grok.shell.addTableView(df);
    await delay(100);
    const renderer = tw.dataFrame.col('SequencesDiff')?.getTag(DG.TAGS.CELL_RENDERER);
    if (renderer !== 'MacromoleculeDifferenceCR')
      throw new Error(`Units 'separator', separator '/' and semType 'MacromoleculeDifference' have been ` +
        `manually set on column but after df aws added as table view renderer has been reset to '${renderer}'`);
  }
});
