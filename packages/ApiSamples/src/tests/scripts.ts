import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, test, delay} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

category('Scripts', () => {
  before(async () => {
    grok.shell.lastError = '';
  });

  const skip = ['function-events', 'demo', 'ui-events'];
  grok.dapi.scripts.filter('package.shortName = "ApiSamples"').list().then((scripts: DG.Script[]) => {
    for (let script of scripts) {
      test(script.friendlyName, async () => {
        // const file = await _package.files.list('scripts', true, script.friendlyName);
        // console.log(file);
        await script.apply();
        await delay(300);
        // if (script.friendlyName === 'function-events') {
        //   grok.functions.onBeforeRunAction.subscribe((_) => {});
        //   grok.functions.onAfterRunAction.subscribe((_) => {});
        // } else if (script.friendlyName === 'demo' || script.friendlyName === 'ui-events') {
        //   grok.events.onEvent('d4-current-view-changed').subscribe((_) => {});
        // }
        if (grok.shell.lastError) {
          const err = grok.shell.lastError;
          grok.shell.lastError = '';
          throw new Error(err);
        }
      }, skip.includes(script.friendlyName) ? {skipReason: 'skip'} : undefined);
    }
  });
});
