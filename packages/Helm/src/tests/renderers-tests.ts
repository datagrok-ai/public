import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {category, expect, test, awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {TAGS as helmTAGS} from '../constants';

category('renderers', () => {
  test('missedInLib', async () => { await _testMissedInLib(); });
});

async function _testMissedInLib() {
  const df: DG.DataFrame = await grok.dapi.files.readCsv('System:AppData/Helm/tests/sample_HELM-missedInLib.csv');
  const helmCol: DG.Column<string> = df.getCol('HELM');

  const tv: DG.TableView = grok.shell.addTableView(df);
  await awaitCheck(() => {
    return $(tv.root).find('.d4-grid canvas').length > 0;
  }, 'Table view canvas not found', 100);

  expect(helmCol.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(helmCol.getTag(DG.TAGS.UNITS), NOTATION.HELM);
  expect(helmCol.getTag(DG.TAGS.CELL_RENDERER), 'helm');

  const cellRendererErrorJson: string = helmCol.getTag(helmTAGS.cellRendererRenderError);
  if (!!cellRendererErrorJson) {
    const cellRendererError: any = JSON.parse(cellRendererErrorJson);
    const err = new Error(cellRendererError['message']);
    err.stack = cellRendererError['stack'];
    throw err;
  }
}
