import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {category, test, testViewer, before, delay} from '@datagrok-libraries/utils/src/test';

category('Viewers', () => {
  const df = grok.data.demo.demog(100);
  const viewers = DG.Func.find({package: 'Charts', tags: ['viewer']}).map((f) => f.friendlyName);
  const viewersToSkip: {[v: string]: string} = {
    'Tree': 'GROK-12569',
    'Word cloud': 'GROK-13198',
    'Sunburst': 'GROK-13778',
  };

  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, await (async () => {
        if (['Sankey', 'Chord'].includes(v)) return (await grok.data.getDemoTable('energy_uk.csv'));
        else if (['Tree', 'Sunburst'].includes(v)) return (await grok.data.getDemoTable('demog.csv'));
        else if (v === 'Globe') return (await grok.data.getDemoTable('geo/earthquakes.csv'));
        return df.clone();
      })(), {detectSemanticTypes: true});
    }, v in viewersToSkip ? {skipReason: viewersToSkip[v]} : {});
  }
});

category('Demo', () => {
  before(async () => {
    grok.shell.lastError = '';
  });

  const demos = DG.Func.find({package: 'Charts', meta: {'demoPath': null}});
  for (const demo of demos) {
    test(demo.friendlyName, async () => {
      await demo.apply();
      await delay(demo.friendlyName === 'sunburstViewerDemo' ? 2000 : 500);
      if (grok.shell.lastError) {
        const err = grok.shell.lastError;
        grok.shell.lastError = '';
        throw new Error(err);
      }
    });
  }
});
