import * as DG from 'datagrok-api/dg';
//import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';

import {category, test, testViewer} from '@datagrok-libraries/utils/src/test';


category('Viewers', () => {
  const viewers = DG.Func.find({package: 'PhyloTreeViewer', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), true);
    }, {skipReason: 'GROK-11534'});
  }
});


const df = DG.DataFrame.fromCsv(
  `node,parent,leaf,distance
root,,false,
node-all,root,false,0.25
node-l1-l2-l3-l4-l5,node-all,false,0.7087542414665222
node-l1-l2-l3,node-l1-l2-l3-l4-l5,false,0.174010768532753
node-l1-l2,node-l1-l2-l3,false,0.17082394659519196
leaf1,node-l1-l2,true,0.13338471949100494
leaf2,node-l1-l2,true,0.051932740956544876
leaf3,node-l1-l2-l3,true,0.18023443222045898
node-l4-l5,node-l1-l2-l3-l4-l5,false,0.20558014512062073
leaf4,node-l4-l5,true,0.4284568130970001
leaf5,node-l4-l5,true,0.6251338720321655
node-l6-l7-l8-l9-l10-l11,node-all,false,0.2287542223930359
node-l6-l7-l8,node-l6-l7-l8-l9-l10-l11,false,0.5740107893943787
node-l6-l7,node-l6-l7-l8,false,0.39653950929641724
leaf6,node-l6-l7,true,0.25338470935821533
leaf7,node-l6-l7,true,0.019327402114868164
leaf8,node-l6-l7-l8,true,0.48023444414138794
node-l9-l10-l11,node-l6-l7-l8-l9-l10-l11,false,0.12234324216842651
node-l9-l10,node-l9-l10-l11,false,0.5755801200866699
leaf9,node-l9-l10,true,0.5684568285942078
leaf10,node-l9-l10,true,0.2451338917016983
leaf11,node-l9-l10-l11,true,0.45654869079589844
leaf-side,root,true,0.3499999940395355`);
