import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function getDemoHitsTargets(): DG.DataFrame {
  let t = grok.shell.table('result');
  let hitsTargets = [];
  for (let i = 0; i < t.rowCount; i++) {
    let targets = Math.floor(Math.random() * 10);
    for (let j = 0; j < targets; j++)
      hitsTargets.push({molregno: t.col('hit')!.get(i), target: 'TRG-' + Math.floor(Math.random() * 100)});
  }
  return DG.DataFrame.fromObjects(hitsTargets)!;
}