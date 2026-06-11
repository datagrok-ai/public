import * as DG from 'datagrok-api/dg';
import {category, test} from '@datagrok-libraries/test/src/test';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';

const FUNC = 'LibTests:TestFileInput';
const FILE_NAME = 'test.txt';
const bytes = () => new TextEncoder().encode('compute-utils file input roundtrip');

// saveRun must accept both compute generations' file-input value types:
// compute1 stores a browser File, compute2 stores a DG.FileInfo.
async function roundtrip(makeValue: () => any) {
  // run twice to ensure repeated save/load works without a DB reset
  for (let i = 0; i < 2; i++) {
    const expected = bytes();
    const fc = DG.Func.byName(FUNC).prepare();
    fc.inputs['inputFile'] = makeValue();

    const saved = await historyUtils.saveRun(fc);
    try {
      const loaded = await historyUtils.loadRun(saved.id);
      const back = await (loaded.inputs['inputFile'] as DG.FileInfo).readAsBytes();
      expectDeepEqual(Array.from(back), Array.from(expected), {prefix: `iteration ${i}`});
    } finally {
      await historyUtils.deleteRun(saved);
    }
  }
}

category('ComputeUtils: History file inputs', async () => {
  test('Save/load file input as browser File (compute1)', async () => {
    await roundtrip(() => new File([bytes()], FILE_NAME));
  });

  test('Save/load file input as DG.FileInfo (compute2)', async () => {
    await roundtrip(() => DG.FileInfo.fromBytes(FILE_NAME, bytes()));
  });
});
