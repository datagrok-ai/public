import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {importFasta, multipleSequenceAlignmentAny} from '../package';

category('renderers', () => {
  let tvList: DG.TableView[];
  let dfList: DG.DataFrame[];

  before(async () => {
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

  async function _testAfterMsa() {
    const fastaTxt: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.fasta');
    const df: DG.DataFrame = importFasta(fastaTxt)[0];
    // await grok.data.detectSemanticTypes(df);

    const srcSeqCol: DG.Column | null = df.col('sequence');
    expect(srcSeqCol !== null, true);
    console.log('Bio: tests/renderers/afterMsa, src data loaded');

    const tv: DG.TableView = grok.shell.addTableView(df);
    console.log('Bio: tests/renderers/afterMsa, table view');

    // await grok.data.detectSemanticTypes(df);
    console.log('Bio: tests/renderers/afterMsa, detectSemanticTypes');

    console.log('Bio: tests/renderers/afterMsa, src before test ' +
      `semType="${srcSeqCol!.semType}", units="${srcSeqCol!.getTag(DG.TAGS.UNITS)}", ` +
      `cell.renderer="${srcSeqCol!.getTag('cell.renderer')}"`);
    expect(srcSeqCol!.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(srcSeqCol!.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
    // TODO: Find the way to check renderer for columns
    expect(srcSeqCol!.getTag('cell.renderer'), 'Macromolecule');
    console.log('Bio: tests/renderers/afterMsa, src semType tested');

    const msaSeqCol: DG.Column | null = await multipleSequenceAlignmentAny(df, srcSeqCol!);
    console.log('Bio: tests/renderers/afterMsa, msaSeqCol created');

    tv.grid.invalidate();
    console.log('Bio: tests/renderers/afterMsa, tv.grid invalidated');

    expect(msaSeqCol!.semType, DG.SEMTYPE.MACROMOLECULE);
    expect(msaSeqCol!.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:PT');
    expect(msaSeqCol!.getTag('cell.renderer'), 'Macromolecule');
    console.log('Bio: tests/renderers/afterMsa, msa semType tested');

    dfList.push(df);
    tvList.push(tv);
  }
});
