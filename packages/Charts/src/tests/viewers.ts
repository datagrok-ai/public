import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {energyUK, demog} from './test-data';


category('Viewers', () => {
  const df = grok.data.demo.demog(100);
  const viewers = DG.Func.find({package: 'Charts', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, (() => {
        if (['SankeyViewer', 'ChordViewer'].includes(v)) return energyUK.clone();
        else if (['TreeMapViewer', 'SunburstViewer'].includes(v)) return demog.clone();
        return df.clone();
      })());
    }, {skipReason: 'GROK-11534'});
  }
});
