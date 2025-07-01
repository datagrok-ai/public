import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, assure, delay, before} from '@datagrok-libraries/utils/src/test';
// import {_package} from '../package-test';
import * as api from '../package-api';
category('Tutorials App', () => {
  test('Launch app', async () => {
    await api.funcs.trackOverview();
    assure.notNull(document.querySelector('div.panel-content > div.tutorials-root'));
  });
});
