import * as DG from 'datagrok-api/dg';
import {Scope} from '@datagrok-libraries/u2';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {Arm, ArmsInput, MAX_ARMS, armsFrame} from '../arms-input';
import {MAX_ARMS as ENGINE_MAX, designStudy, tooManyArms} from '../engine';
import {declarationOf} from '../frames';
import {CONTROL_COLOR, armColors} from '../colors';
import {Band, MANY_BANDS, MAX_LEGEND_LINES, chartOption, legendLayout, legendLines} from '../band-chart';

const arms = (k: number): Arm[] => Array.from({length: k}, (_, i) =>
  ({name: `Cpd-${i + 1}`, therapy: 'small_molecule', metric: 'TGI', effect: 50 + 2 * i}));
const base = {line: 'MC38', endpoint: 1000, power: 80, correction: 'holm'};

// Up to 15 treated arms (owner, 0.1.3): planned whole, 16 colours apart, the chart readable, the 16th refused.
category('Fifteen arms', () => {
  test('the limit is one: 15 in the engine and the arms input', async () => {
    expect(MAX_ARMS, 15);
    expect(ENGINE_MAX, MAX_ARMS);
  });

  test('a 16th arm is refused with the reason, before any engine call', async () => {
    let said = '';
    try {
      await designStudy(declarationOf({...base, arms: armsFrame(arms(16))} as any));
    } catch (e: any) {
      said = String(e?.message ?? e);
    }
    expect(said, tooManyArms(16));
    expect(said.includes('at most 15'), true, said);
  });

  test('the arms input stops at 15: Add treatment arm says why', async () => {
    const root = Scope.root('fifteen arms test');
    try {
      const input = root.build(() => new ArmsInput({}));
      input.value = armsFrame(arms(15));
      await delay(50);
      input.add();
      await delay(50);
      expect(input.count.value, 15);
      const add = input.root.querySelector<HTMLButtonElement>('[data-u2=pmx-add-arm]')!;
      expect(add.title, 'At most 15 treated arms');
    } finally {
      root.dispose();
    }
  });

  test('15 treated arms and the control: 16 colours, none twice, no treated arm grey', async () => {
    const colors = armColors(arms(15).map((a) => a.name));
    const all = Object.values(colors);
    expect(all.length, 16);
    expect(new Set(all.map((c) => c.toLowerCase())).size, 16, all.join(' '));
    const grey = (hex: string) => {
      const c = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16));
      return Math.max(...c) - Math.min(...c) < 24;
    };
    expect(all.filter((c) => c !== CONTROL_COLOR).some(grey), false, all.join(' '));
  });

  test('the growth chart past 7 groups: one range at a time, a legend wrapped, never dropped', async () => {
    const band = (name: string): Band => ({name, x: [0, 1], lower: [1, 2], upper: [3, 4], line: [2, 3], row: 0});
    const bands = ['Control', ...arms(15).map((a) => a.name)].map(band);
    const look = (focus: string | null) => ({xTitle: 'Day', yTitle: 'Volume', log: false, markers: false, fit: false,
      dashed: [], legend: false, xPad: false, focus, fmt: {xUnits: '', yUnits: '', x: '0', y: '0', xTitle: 'Day'}});
    const drawn = (focus: string | null) => (chartOption(bands, bands.map(() => '#1f77b4'), [], look(focus))
      .series as any[]).filter((s) => s.type === 'custom')
      .map((s) => s.renderItem({dataIndex: 0}, {coord: (p: number[]) => p}) !== null);
    expect(drawn(null).filter((x) => x).length, 0, 'none pointed at: the lines alone');
    expect(drawn('Cpd-3').map((x, i) => x ? bands[i].name : '').join(''), 'Cpd-3');
    const few = bands.slice(0, MANY_BANDS);
    expect((chartOption(few, few.map(() => '#1f77b4'), [], look(null)).series as any[])
      .filter((s) => s.type === 'custom')
      .every((s) => s.renderItem({dataIndex: 0}, {coord: (p: number[]) => p}) !== null), true, '7 groups: every range');
    const names = bands.map((b) => b.name);
    expect(legendLines(names, 2000), 1);
    expect(legendLines(names, 600) > 1, true);
    const wrapped = chartOption(bands, bands.map(() => '#1f77b4'), [], {...look(null), legend: true,
      legendLayout: legendLayout(names, 1100)});
    expect((wrapped.legend as any).show, true, '16 groups: a legend still');
    expect((wrapped.legend as any).type, 'plain', 'wrapped, not paged');
    expect((wrapped.grid as any).top, 8 + 24 * legendLines(names, 1100), 'the plot below its lines');
    expect(legendLines(names, 1100) <= MAX_LEGEND_LINES, true, String(legendLines(names, 1100)));
    // past a few lines (12 studies in Compare): one scrolled line, the plot never pushed away (chaos round 3 #22)
    expect(JSON.stringify(legendLayout(names, 200)), JSON.stringify({lines: 1, type: 'scroll'}));
  });

  test('15 treated arms planned end to end: every arm sized, in the tables and the curves', async () => {
    const fc = DG.Func.find({package: 'PmxStudio', name: 'pmxPlanStudy'})[0].prepare({...base,
      arms: armsFrame(arms(15))});
    await fc.call();
    const armResults = fc.outputs['armResults'] as DG.DataFrame;
    const plan = fc.outputs['plan'] as DG.DataFrame;
    const curves = fc.outputs['curves'] as DG.DataFrame;
    expect(armResults.rowCount, 15);
    const n = Array.from({length: 15}, (_, i) => armResults.get('n', i) as number);
    expect(n.every((x) => Number.isFinite(x) && x > 0), true, n.join(','));
    expect(plan.get('n_control', 0), Math.max(...n));
    const groups = new Set(curves.col('arm')!.toList());
    expect(groups.size, 16, [...groups].join(','));
  }, {timeout: 900000});
});
