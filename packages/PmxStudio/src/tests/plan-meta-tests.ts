import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SavedPlan} from '../plan-file';
import {StudyHolds, planOnStudy} from '../plan-vs-observed';
import {StudyMeta, emptyMeta} from '../study-meta';

// a plan with two arms and its study's details, as Plan saves them
const plan = (study: Partial<StudyMeta> = {}, inputs: {[k: string]: unknown} = {}): SavedPlan => ({
  format: 'pmx-plan', version: 1, name: 'CRC-021', created: '', author: '',
  study: {...emptyMeta(), runBy: 'Ana', start: '2026-12-01', colors: {'Drug A 10 mg/kg': '#112233'},
  groups: [{group: 'Drug A 10 mg/kg', compound: 'Drug A', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'QD'},
    {group: 'Drug B 5 mg/kg', compound: 'Drug B', dose: 5, unit: 'mg/kg', route: 'p.o.', schedule: 'BID'}],
  ...study},
  inputs: {looks: '7, 14', measurementDays: 'Mon, Thu', arms: [{name: 'Drug A 10 mg/kg'}, {name: 'Drug B 5 mg/kg'}],
    ...inputs},
  predicted: {line: 'MC38', endpoint: 1500},
} as unknown as SavedPlan);

const holds = (x: Partial<StudyHolds> = {}): StudyHolds => ({meta: null, looks: '', line: null, plan: null,
  groups: ['Vehicle', 'drug a 10 MG/KG', 'Drug C'], ...x});

category('Plan picked: the study\'s details', () => {
  test('fills what the study lacks, groups matched by name in any case', async () => {
    const r = planOnStudy(plan(), holds());
    expect(r.differences.length, 0);
    expect(r.looks, '7, 14');
    expect(r.line, 'MC38');
    expect(r.carry, true);
    expect(r.meta?.runBy, 'Ana');
    expect(r.meta?.start, '2026-12-01');
    // the plan's row goes under the data's own name for the group; Drug B (no group) is not invented
    expect(r.meta?.groups.length, 1);
    expect(r.meta?.groups[0].group, 'drug a 10 MG/KG');
    expect(r.meta?.groups[0].schedule, 'QD');
    expect(r.meta?.colors['drug a 10 MG/KG'], '#112233');
    expect(r.missing.join(), 'Drug B 5 mg/kg');
    expect(r.extra.join(), 'Drug C');
  });

  test('a value the study holds otherwise is kept, and said', async () => {
    const meta = {...emptyMeta(), runBy: 'Ana', start: '2026-11-30', groups: [{group: 'drug a 10 MG/KG',
      compound: 'Drug A', dose: 20, unit: 'mg/kg', route: 'i.p.', schedule: 'QD'}]};
    const s = holds({meta, looks: '7,14', line: 'CT26'});
    const r = planOnStudy(plan(), s);
    // the reviews read the same (7, 14 written otherwise); Run by is the same
    // and the arms that do not match, with how many (QA 0.1.3 round 4 R4-2: the model alone was listed)
    expect(r.differences.length, 4);
    expect(r.differences.includes('Treated groups: the plan has 2 (Drug A 10 mg/kg, Drug B 5 mg/kg), the study ' +
      'has 2 (drug a 10 MG/KG, Drug C)'), true, r.differences.join(' | '));
    expect(r.differences.some((d) => d.startsWith('Day 0: the plan says 2026-12-01, the study has 2026-11-30')),
      true);
    expect(r.differences.some((d) => d.startsWith('Tumor model:')), true);
    expect(r.differences.some((d) => d.startsWith('drug a 10 MG/KG: the plan says Drug A 10 mg/kg')), true);
    expect(r.meta?.start, '2026-11-30');
    expect(r.meta?.groups[0].dose, 20);
    expect(r.line, 'CT26');
    expect(r.looks, '7,14');
    // with the plan's taken
    const u = planOnStudy(plan(), s, true);
    expect(u.meta?.start, '2026-12-01');
    expect(u.meta?.groups[0].dose, 10);
    expect(u.meta?.groups.length, 1);
    expect(u.line, 'MC38');
  });

  test('a study carrying another plan: its measurement days are asked about', async () => {
    const other = plan({}, {measurementDays: 'Tue, Fri'});
    other.name = 'Another';
    const r = planOnStudy(plan(), holds({plan: other}));
    expect(r.carry, false);
    expect(r.differences.includes('Measurement days: the plan says Mon, Thu, the study has Tue, Fri'), true);
    expect(planOnStudy(plan(), holds({plan: other}), true).carry, true);
    // the same plan: nothing to ask, nothing carried again
    const same = planOnStudy(plan(), holds({plan: plan()}));
    expect(same.carry, false);
    expect(same.differences.length, 0);
  });

  test('a plan without study details leaves the study\'s as they are', async () => {
    const p = plan();
    delete p.study;
    expect(planOnStudy(p, holds()).meta, null);
    const meta = {...emptyMeta(), runBy: 'Bo'};
    expect(planOnStudy(p, holds({meta})).meta, meta);
  });
});

category('Plan picked: reviews fixed', () => {
  test('fixed reviews are neither asked about nor changed', async () => {
    const r = planOnStudy(plan(), holds({looks: '10', looksFixed: true}));
    expect(r.looks, '10');
    expect(r.differences.length, 0);
  });
});
