import {before, category, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from './utils';

category('UI: List', () => {
  let v: DG.View;

  const list = ui.list([
    'element 1',
    'element 2',
  ]);

  before(async () => {
    v = grok.shell.newView('');
  });

  test('list.root', async () => {
    checkHTMLElement('list', list, v, '.d4-flex-col');
  });
});
