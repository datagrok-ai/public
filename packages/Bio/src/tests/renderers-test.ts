import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {importFasta, multipleSequenceAlignmentAny} from '../package';
import {readDataframe} from './utils';
import {convertDo} from '../utils/convert';
import {ALPHABET, NOTATION, UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

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

  async function _testAfterMsa() {
    const fastaTxt: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.fasta');
    const df: DG.DataFrame = importFasta(fastaTxt)[0];
    const tv: DG.TableView = grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);
    console.log('Bio: tests/renderers/afterMsa, table view');

    const srcSeqCol: DG.Column | null = df.col('sequence');
    expect(srcSeqCol !== null, true);

    console.log('Bio: tests/renderers/afterMsa, src before test ' +
      `semType="${srcSeqCol!.semType}", units="${srcSeqCol!.getTag(DG.TAGS.UNITS)}", ` +
      `cell.renderer="${srcSeqCol!.getTag(DG.TAGS.CELL_RENDERER)}"`);
    expect(srcSeqCol!.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(srcSeqCol!.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(srcSeqCol!.getTag(UnitsHandler.TAGS.aligned), 'SEQ');
    expect(srcSeqCol!.getTag(UnitsHandler.TAGS.alphabet), ALPHABET.PT);
    expect(srcSeqCol!.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    const msaSeqCol: DG.Column | null = await multipleSequenceAlignmentAny(df, srcSeqCol!);
    tv.grid.invalidate();

    expect(msaSeqCol!.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(msaSeqCol!.getTag(DG.TAGS.UNITS), NOTATION.FASTA);
    expect(msaSeqCol!.getTag(UnitsHandler.TAGS.aligned), 'SEQ.MSA');
    expect(msaSeqCol!.getTag(UnitsHandler.TAGS.alphabet), ALPHABET.PT);
    expect(msaSeqCol!.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    dfList.push(df);
    tvList.push(tv);
  }

  async function _testAfterConvert() {
    const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA_PT.csv');
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    const tv: DG.TableView = grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);

    const srcCol: DG.Column = df.col('sequence')!;
    const tgtCol: DG.Column = await convertDo(srcCol, NOTATION.SEPARATOR, '/');
    expect(tgtCol.getTag(DG.TAGS.CELL_RENDERER), 'sequence');

    tvList.push(tv);
    dfList.push(df);
  };
});
