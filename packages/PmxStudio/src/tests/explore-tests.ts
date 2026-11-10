import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {AnalysisFrames, exampleStudy} from '../analyze';
import {armColors} from '../frames';
import {COL, RESPONSE, StudyInputs, animalsTable, explore} from '../explore';

const OUTPUTS = ['analysis', 'comparisons', 'animals', 'volumes', 'twins', 'effects', 'checks', 'versusPlan',
  'predicted', 'tolerability', 'weights'];

async function example(): Promise<{r: AnalysisFrames, s: StudyInputs, colors: {[g: string]: string}}> {
  const study = await exampleStudy();
  const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxAnalyzeStudy'})[0].prepare({study,
    animal: study.col('animal'), group: study.col('group'), day: study.col('day'), volume: study.col('volume'),
    line: 'MC38', endpoint: 1500, correction: 'holm'});
  await fc.call();
  const r = Object.fromEntries(OUTPUTS.map((o) => [o, fc.outputs[o]])) as unknown as AnalysisFrames;
  const colors = armColors(r.comparisons.col('arm')!.toList());
  colors[r.analysis.get('control', 0)] = colors['Control'];
  return {r, colors, s: {study, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null}};
}

category('Your data', () => {
  test('the animals: one row each, the response as the twin reads it', async () => {
    const {r, s, colors} = await example();
    const a = animalsTable(r, s, colors);
    expect(a.rowCount, r.animals.rowCount);
    const responders = (a.col(COL.response)!.toList() as string[]).filter((x) => x === RESPONSE.responder).length;
    expect(responders, (r.effects.col('responder')!.toList() as boolean[]).filter((x) => x).length);
    const controls = (a.col(COL.group)!.toList() as string[]).filter((g) => g === r.analysis.get('control', 0));
    expect(controls.length, r.analysis.get('n_control', 0));
    // every animal has its volume at randomization
    expect(a.col(COL.v0)!.stats.missingValueCount, 0);
  }, {timeout: 120000});

  test('the measurements: the study\'s rows, linked to the animals by selection', async () => {
    const {r, s, colors} = await example();
    explore('animals', r, s, colors, grok.shell.v);
    await new Promise((res) => setTimeout(res, 50));
    const animals = grok.shell.tv.dataFrame;
    expect(animals.name.startsWith('Animals'), true, animals.name);
    const types = [...grok.shell.tv.viewers].map((v) => v.type);
    expect(types.includes(DG.VIEWER.BOX_PLOT) && types.includes(DG.VIEWER.SCATTER_PLOT), true, types.join());
    const measurements = [...grok.shell.tableViews].map((v) => v.dataFrame)
      .find((df) => df.name.startsWith('Measurements'))!;
    expect(measurements.rowCount, s.study.rowCount);
    // one animal selected: its measurements selected
    const id = animals.get(COL.animal, 0);
    animals.selection.init((i) => i === 0);
    await new Promise((res) => setTimeout(res, 200));
    const own = (measurements.col(COL.animal)!.toList() as string[]).filter((x) => x === id).length;
    expect(measurements.selection.trueCount, own);
    for (const v of [...grok.shell.tableViews].filter((x) => [animals, measurements].includes(x.dataFrame)))
      v.close();
  }, {timeout: 120000});
});
