import * as grok from 'datagrok-api/grok';

import {before, category, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {getBaseURL, setPackage} from '../package';
import * as api from '../package-api';

category('Alation', () => {
  before(async () => setPackage(_package));

  test('App test', async () => {
    try {
      await getBaseURL();
    } catch (_) {
      return;
    }
    await api.funcs.alation();
  }, {skipReason: 'Running tests manually'});
});
