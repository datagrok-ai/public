import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {CONTROL_COLOR, StudyLink, armColors, follow, groupColors, hasOwnColors, linkStudy, paint, recolorOf,
  resetColors, slotsAfter, treatedColors} from '../colors';
import {StudyMeta, emptyMeta, f32, parseMeta, renamedDetails} from '../study-meta';
import {openStudyTable, tagStudy} from '../study-table';

const same = (a: unknown, b: unknown, what: string) => expect(JSON.stringify(a), JSON.stringify(b), what);

/** A stage open on a study: its details kept in memory, its groups `groups` with `control`. */
function stage(id: string, groups: string[], control: string, colors: {[g: string]: string} = {}):
  StudyLink & {kept: StudyMeta} {
  const s = {kept: {...emptyMeta(), id, colors}} as StudyLink & {kept: StudyMeta};
  s.meta = () => s.kept;
  s.setMeta = (m) => s.kept = m;
  s.colors = () => groupColors(groups, control, s.kept.colors);
  return s;
}

const groupFrame = (groups: string[]) => DG.DataFrame.fromColumns([DG.Column.fromStrings('group', groups)]);
const colorOf = (df: DG.DataFrame, g: string) =>
  JSON.parse(df.col('group')!.getTag(DG.TAGS.COLOR_CODING_CATEGORICAL) || '{}')[g];

category('Study colours', () => {
  test('an arm keeps its colour for life: a combination added, one renamed or removed recolours no other', async () => {
    const [a, b, ab] = ['LD4172 20 mg/kg', 'anti-PD-1 10 mg/kg', 'LD4172 20 mg/kg + anti-PD-1 10 mg/kg'];
    const two = slotsAfter({}, [], [a, b]);
    const three = slotsAfter(two, [a, b], [a, b, ab]);
    const before = treatedColors([a, b], two);
    const after = treatedColors([a, b, ab], three);
    same([after[a], after[b]], [before[a], before[b]], 'the first two keep theirs');
    expect(new Set(Object.values(after)).size, 3, 'the combination takes a free one');
    // Monitor and Analyze list the groups in another order: the same colours
    same(groupColors(['Vehicle', ab, b, a], 'Vehicle', {}, three)[a], after[a], 'in any stage');
    const renamed = slotsAfter(three, [a, b, ab], [a, 'anti-PD-1 5 mg/kg', ab]);
    expect(renamed['anti-PD-1 5 mg/kg'], three[b], 'renamed where it stands');
    same(slotsAfter(three, [a, b, ab], [a, ab]), {[a]: three[a], [ab]: three[ab]}, 'one removed');
    const sorted = (c: {[g: string]: string}) => Object.entries(c).sort();
    same(sorted(treatedColors([b, a])), sorted(treatedColors([a, b])), 'no places yet: by name, whatever the order');
  });

  test('a dose read back from a 32-bit column is the dose typed (0.3, not 0.30000001…)', async () => {
    const c = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'dose', [0.3, 14.96, 100]);
    expect(f32(c.get(0)) === 0.3, true);
    expect(f32(c.get(1)) === 14.96, true);
    expect(f32(c.get(2)) === 100, true);
  });
  test('a group without a colour of its own has the automatic one; one with, its own', async () => {
    const auto = armColors(['B', 'A']);
    same(auto, {'Control': CONTROL_COLOR, ...treatedColors(['B', 'A'])}, 'the automatic palette');
    const own = armColors(['B', 'A'], {'A': '#112233', 'Gone': '#445566'});
    expect(own['A'], '#112233', 'the user\'s colour');
    expect(own['B'], auto['B'], 'the others automatic');
    expect('Gone' in own, false, 'a colour of a group not in the study is not listed');
    // a new group: the automatic colour of its name
    expect(armColors(['B', 'A', 'C'], {'A': '#112233'})['C'], treatedColors(['B', 'A', 'C'])['C'], 'a new group');
    const g = groupColors(['Vehicle', 'A'], 'Vehicle', {'Vehicle': '#000000'});
    expect(g['Vehicle'], '#000000', 'a control of its own colour');
    expect(groupColors(['Vehicle', 'A'], 'Vehicle')['Vehicle'], CONTROL_COLOR, 'the control in its slate');
  });

  test('old files have no colours: the automatic palette; a bad colour is dropped', async () => {
    same(parseMeta({id: 'PMX-1'})!.colors, {}, 'no colours');
    same(parseMeta({colors: {'A': '#aabbcc', 'B': 'red', 'C': 3}})!.colors, {'A': '#aabbcc'}, 'only hex colours');
    same(parseMeta(JSON.parse(JSON.stringify({...emptyMeta(), colors: {'A': '#aabbcc'}})))!.colors,
      {'A': '#aabbcc'}, 'a round trip through the file');
  });

  test('a renamed group keeps its colour', async () => {
    const m = {...emptyMeta(), colors: {'A': '#aabbcc', 'B': '#000000'}};
    same(renamedDetails(m, new Map([['A', 'A2']])).colors, {'B': '#000000', 'A2': '#aabbcc'}, 'moved');
    same(renamedDetails(m, new Map([['A', 'B']])).colors, {'B': '#000000'}, 'the new name\'s own colour wins');
  });

  test('a colour changed in one chart is the study\'s, and its other stages follow', async () => {
    const groups = ['Vehicle', 'A', 'B'];
    const plan = stage('PMX-T-COLOR', groups, 'Vehicle');
    const monitor = stage('PMX-T-COLOR', groups, 'Vehicle');
    const other = stage('PMX-T-OTHER', ['Vehicle', 'X'], 'Vehicle');
    const stops = [plan, monitor, other].map(linkStudy);
    try {
      const shown = groupFrame(groups);
      paint(shown, 'group', plan.colors());
      const elsewhere = groupFrame(groups);
      follow(monitor, [elsewhere]);
      expect(colorOf(elsewhere, 'A') === undefined, true, 'a column not coloured by group is left as it is');
      paint(elsewhere, 'group', monitor.colors());
      await delay(10);
      expect(hasOwnColors(plan), false, 'no colour of its own yet: our own painting is not the user\'s');
      // the user picks another colour for A in a viewer of the first frame
      shown.col('group')!.meta.colors.setCategorical({...plan.colors(), 'A': '#123456'});
      await delay(100);
      expect(plan.kept.colors['A'], '#123456', 'the study keeps it');
      expect(monitor.kept.colors['A'], '#123456', 'the same study in another stage');
      expect(other.kept.colors['A'] === undefined, true, 'another study does not');
      follow(monitor, [elsewhere]);
      expect(colorOf(elsewhere, 'A'), '#123456', 'the other stage\'s frame follows');
      expect(colorOf(elsewhere, 'B'), monitor.colors()['B'], 'the other groups as they were');
      resetColors(plan);
      same([plan.kept.colors, monitor.kept.colors], [{}, {}], 'reset: the automatic colours, in every stage');
      follow(monitor, [elsewhere]);
      expect(colorOf(elsewhere, 'A'), treatedColors(['A', 'B'])['A'], 'and its frames');
    } finally {
      stops.forEach((s) => s());
    }
  });

  test('Plan, Monitor and Analyze of one study share its colours, the control\'s too', async () => {
    const PLAN = '["MC38 full study",{"arms":[]}]';
    // a study whose plan has no ID yet: Plan knows it by its plan; Plan's control arm is 'Control', the study's
    // control group 'Vehicle'
    const at = (groups: string[], control: string) => {
      const s = stage('', groups, control);
      s.plan = () => PLAN;
      s.control = () => control;
      return s;
    };
    const monitor = at(['Vehicle', 'anti-PD-1'], 'Vehicle');
    const analyze = at(['Vehicle', 'anti-PD-1'], 'Vehicle');
    const plan = at(['Control', 'anti-PD-1'], 'Control');
    const stops = [monitor, analyze, plan].map(linkStudy);
    let late: (() => void) | null = null;
    try {
      recolorOf('anti-PD-1', monitor.colors()['anti-PD-1'])!('#ff0000');
      expect(plan.colors()['anti-PD-1'], '#ff0000', 'Plan, open before the change, follows Monitor');
      expect(analyze.colors()['anti-PD-1'], '#ff0000', 'Analyze too');
      recolorOf('Vehicle', monitor.colors()['Vehicle'])!('#00ff00');
      expect(plan.colors()['Control'], '#00ff00', 'the control recoloured in Monitor: Plan\'s control arm');
      expect(analyze.colors()['Vehicle'], '#00ff00', 'and Analyze\'s control');
      recolorOf('Control', plan.colors()['Control'])!('#0000ff');
      expect(monitor.colors()['Vehicle'], '#0000ff', 'the control recoloured in Plan: the study\'s control group');
      const plan2 = at(['Control', 'anti-PD-1'], 'Control');
      late = linkStudy(plan2);
      expect(plan2.colors()['anti-PD-1'], '#ff0000', 'Plan opened after the change takes it');
      expect(plan2.colors()['Control'], '#0000ff', 'and the control\'s');
    } finally {
      stops.forEach((s) => s());
      late?.();
    }
  });

  test('two open studies with a group of one name and colour: the colour picked is the table\'s study\'s', async () => {
    const [one, two] = ['PMX-T-ONE', 'PMX-T-TWO'].map((id) => stage(id, ['Vehicle', 'A'], 'Vehicle'));
    const [t1, t2] = [groupFrame(['A']), groupFrame(['A'])];
    one.table = () => t1;
    two.table = () => t2;
    const stops = [linkStudy(one), linkStudy(two)];
    try {
      expect(one.colors()['A'], two.colors()['A'], 'the same colour in both');
      recolorOf('A', two.colors()['A'], t2)!('#123123');
      expect([one.kept.colors['A'] ?? 'none', two.kept.colors['A']].join(), 'none,#123123');
    } finally {
      stops.forEach((s) => s());
    }
  });

  test('a colour picked on a chart\'s legend item is the study\'s', async () => {
    const s = stage('PMX-T-LEGEND', ['Vehicle', 'A'], 'Vehicle');
    const stop = linkStudy(s);
    try {
      expect(recolorOf('A', '#000001') === null, true, 'no study draws A in that colour');
      recolorOf('A', s.colors()['A'])!('#ABCDEF');
      expect(s.kept.colors['A'], '#abcdef');
    } finally {
      stop();
    }
  });

  test('the study\'s table opens in the study\'s colours', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['a', 'b']),
      DG.Column.fromStrings('group', ['Vehicle', 'A']), DG.Column.fromList(DG.TYPE.FLOAT, 'day', [0, 0]),
      DG.Column.fromList(DG.TYPE.FLOAT, 'volume', [100, 110])]);
    tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null});
    const s = stage('PMX-T-TABLE', ['Vehicle', 'A'], 'Vehicle', {'A': '#123456'});
    s.table = () => df;
    const stop = linkStudy(s);
    try {
      openStudyTable(df);
      expect(colorOf(df, 'A'), '#123456');
      expect(colorOf(df, 'Vehicle'), CONTROL_COLOR);
    } finally {
      stop();
      [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id)?.close();
    }
  });
});
