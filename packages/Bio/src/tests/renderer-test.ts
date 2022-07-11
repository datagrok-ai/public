import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {importFasta, multipleSequenceAlignmentAny} from '../package';

category('renderers', () => {
  test('afterMsa', async () => {
    await _testAfterMsa();
  });
});

export async function _testAfterMsa() {
  const fastaTxt: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.fasta');
  const df: DG.DataFrame = importFasta(fastaTxt)[0];

  const seqCol: DG.Column | null = df.col('sequence');
  expect(seqCol !== null, true);

  const tv: DG.TableView = grok.shell.addTableView(df);
  await grok.data.detectSemanticTypes(df);

  expect(seqCol!.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(seqCol!.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
  expect(seqCol!.getTag('cell.renderer'), 'Macromolecule');

  const seqMsaCol: DG.Column = await multipleSequenceAlignmentAny(df, seqCol!);
  tv.grid.invalidate();

  expect(seqMsaCol!.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(seqMsaCol!.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:PT');
  expect(seqMsaCol!.getTag('cell.renderer'), 'Macromolecule');

  // tv.close();
}
