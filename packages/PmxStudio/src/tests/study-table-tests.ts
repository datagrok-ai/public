import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_RULES, findOutliers, StudyColumns} from '../outlier-rules';
import {FLAG, FLAG_REASON, appendPasted, emptyStudy, exclusionsOf, withoutExcluded} from '../study-table';

/** Six animals of one group on days 0, 3 and 7; M4's day 7 a caliper typo (2400 for ~240), M1 losing weight. */
function study(): StudyColumns {
  const v: {[a: string]: number[]} = {M1: [100, 150, 220], M2: [105, 155, 230], M3: [95, 145, 210],
    M4: [110, 160, 2400], M5: [98, 148, 215], M6: [102, 152, 225]};
  const s: StudyColumns = {animal: [], group: [], day: [], volume: [], weight: []};
  for (const [a, vs] of Object.entries(v)) {
    vs.forEach((x, k) => {
      s.animal.push(a);
      s.group.push('Vehicle');
      s.day.push([0, 3, 7][k]);
      s.volume.push(x);
      s.weight!.push(a === 'M1' ? [20, 20, 17][k] : 20);
    });
  }
  return s;
}

category('Study table', () => {
  test('the rules find the typo and the weight loss, nothing else', async () => {
    const s = study();
    const found = findOutliers(s, DEFAULT_RULES);
    const typo = s.animal.findIndex((a, i) => a === 'M4' && s.day[i] === 7);
    const thin = s.animal.findIndex((a, i) => a === 'M1' && s.day[i] === 7);
    expect(found.filter((f) => f.row === typo && f.on === 'volume').length, 2, 'jump and far from its group');
    expect(found.filter((f) => f.row === thin && f.on === 'weight').length, 1, '15% lighter');
    expect(found.every((f) => f.row === typo || f.row === thin), true, found.map((f) => f.reason).join('; '));
  });

  test('rules switched off find nothing; impossible values are found', async () => {
    const off = {impossible: false, jump: {on: false, fold: 3}, group: {on: false, mads: 3.5, minAnimals: 5},
      weightDrop: {on: false, percent: 10}};
    expect(findOutliers(study(), off).length, 0);
    const s = study();
    s.volume[0] = -5;
    s.day[1] = 0; // M1 measured twice on day 0
    const reasons = findOutliers(s, {...off, impossible: true}).map((f) => f.reason).join('; ');
    expect(reasons.includes('negative volume'), true, reasons);
    expect(reasons.includes('measured twice'), true, reasons);
  });

  test('pasted rows land in their columns, with or without a header', async () => {
    const df = emptyStudy();
    expect(appendPasted(df, 'Animal\tGroup\tDay\tVolume\r\nM1\tVehicle\t0\t101,5\nM2\tVehicle\t0\t99\n'), 2);
    expect(df.get('volume', 0), 101.5);
    expect(df.get('animal', 1), 'M2');
    expect(df.col('body_weight')!.isNone(0), true);
    expect(appendPasted(df, 'M3\tVehicle\t3\t150\t21'), 1);
    expect(df.get('body_weight', 2), 21);
    expect(df.get('day', 2), 3);
  });

  test('excluded rows leave the engine’s inputs and are reported with their reasons', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M2', 'M2']),
      DG.Column.fromStrings('group', ['A', 'A', 'A', 'A']), DG.Column.fromList('double', 'day', [0, 3, 0, 3]),
      DG.Column.fromList('double', 'volume', [100, 150, 90, 900]),
      DG.Column.fromStrings(FLAG, ['excluded', 'excluded', '', 'excluded']),
      DG.Column.fromStrings(FLAG_REASON, ['dosing error', 'dosing error', '', 'misread'])]);
    const x = withoutExcluded({study: df, animal: df.col('animal')!, endpoint: 1500});
    expect(x.study.rowCount, 1);
    expect(x.animal.length, 1, 'columns come from the kept table');
    expect(x.endpoint, 1500);
    expect(df.rowCount, 4, 'the study itself keeps every row');
    const text = exclusionsOf({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume',
      bodyWeight: null})!;
    expect(text.includes('animal M1, all its measurements (dosing error)'), true, text);
    expect(text.includes('M2 on day 3 (misread)'), true, text);
    const clean = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1'])]);
    expect(withoutExcluded({study: clean}).study === clean, true, 'nothing excluded: the same table');
  });
});
