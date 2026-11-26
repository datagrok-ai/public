import {category, expect, test} from '@datagrok-libraries/test/src/test';

import * as DG from 'datagrok-api/dg';

import {slopeSide} from '../analyze';
import {chartRows, compareFooter, compareRows, controlCheck, controlGrowth, controlVsAtlas, doseGroups, driftSentence,
  growthFrame, rangeText, sideText, compareSection, compareSentence, compareTexts, compareView, hintRows,
  humansToOpen, modelWarning, studyLabels, tgiFootnote, tgiNote, withHumans}
  from '../compare-view';
import {entryOf} from '../library';
import {EXPOSURE_TAG} from '../study-file';
import {emptyMeta, setMetaOn} from '../study-meta';
import {shownMultiple} from '../exposure';
import {Exposure, HumanExposure, NcaByArm, ReferenceRow, exposureJson, exposureRange, fromReference, patientsRange}
  from '../exposure';
import {ExposureSection, GroupDetails, detailsElement, exposureSection, groupDetails} from '../translation-report';
import type {PreviousStudy} from '../study-effect';
import {drugKey, plainArm} from '../arm-name';
import {animalPanel, animalReadings} from '../study-table';
import {sig} from '../sheet';

const x = (dose: number, aucInf: number): Exposure => ({auc: aucInf * 0.9, aucCv: null, aucUnit: 'mg/L·h',
  cmax: null, cmaxCv: null, cmaxUnit: 'mg/L', aucInf, aucInfUnit: 'mg/L·h', dose: `${dose} mg/kg`, route: 'PO',
  schedule: 'QD', profiles: 6, drug: 'Cpd-A', pooled: false});
const patients = (drug: string, h: Partial<HumanExposure> = {}): [string, HumanExposure] =>
  [drugKey(drug), {drug, auc: 20000, cmax: null, unit: 'ng/mL', fuMouse: null, fuHuman: null, source: 'projected',
    cv: 45, ...h}];
/** A study of `groups` (3 animals each on day 14, the control at 1000 mm³), with what it kept. */
const study = (name: string, groups: [string, number][], h: [string, HumanExposure][] = [],
  nca: NcaByArm | null = null): PreviousStudy => ({name, line: null, meta: null,
  rows: [['Vehicle', 1000] as [string, number], ...groups].flatMap(([g, v]) => [1, 2, 3].map((i) =>
    ({animal: `${g}-${i}`, group: g, day: 14, volume: v}))),
  exposure: h.length || nca ? exposureJson(new Map(h), nca) : undefined});

category('Compare: studies against the patients\' exposure', () => {
  const a = study('CRC-017', [['Cpd-A 10 mg/kg QD', 600], ['Cpd-A 30 mg/kg QD', 400]], [patients('Cpd-A')],
    {table: 'pk', pk: new Map(), arms: new Map([['Cpd-A 10 mg/kg QD', {exposure: x(10, 10)}],
      ['Cpd-A 30 mg/kg QD', {exposure: x(30, 40)}]])});
  const b = study('CRC-021', [['Cpd-B 25 mg/kg', 700], ['Cpd-B 75 mg/kg', 500]]);
  // a reference drug without mouse PK: in the table by dose, not on the chart
  const c = study('CRC-019', [['Paclitaxel 20 mg/kg', 450]], [patients('Paclitaxel', {source: 'measured'})]);
  const rows = [...compareRows(a, 'HCT116'), ...compareRows(b, 'HCT116'), ...compareRows(c, 'HCT116')];

  test('each treated group with its TGI against its own control', async () => {
    expect(rows.length, 5);
    expect(rows.map((r) => Math.round(r.tgi!)).join(','), '40,60,30,50,55');
    expect(rows[0].source, 'Projected');
    expect(tgiNote(rows[0]), '1 − T/C of medians, day 14 (not from the analysis)');
    // the analysis's TGI and day, when the study kept them
    const kept = {...a, exposure: {...a.exposure, effect: {day: 21, tgi: {'Cpd-A 10 mg/kg QD': 44.5}}}};
    const [first, second] = compareRows(kept, 'HCT116');
    expect([first.tgi, first.day, first.analyzed, tgiNote(first)].join('|'), '44.5|21|true|');
    expect(second.analyzed, false, 'a group the analysis did not keep: the medians');
    // the medians starred and said so once; none starred: no footnote
    expect(tgiFootnote(compareRows(kept, 'HCT116'))!.startsWith('* The study kept no analysis'), true);
    expect(tgiFootnote([first]), null);
  });

  test('the chart: only groups with a TGI and a multiple of the patients\' exposure', async () => {
    expect(chartRows(rows).map((r) => r.group).join('|'), 'Cpd-A 10 mg/kg QD|Cpd-A 30 mg/kg QD');
    const pac = rows.find((r) => r.drug === 'Paclitaxel')!;
    expect(pac.multiple?.ratio === undefined, true);
    expect(!!pac.multiple?.why, true, 'the reason it is not placed');
    // one range for the band and the arrows: CV 50 % from 0.41 to 1.95, no CV 1/3 to 3
    const band = patientsRange(null);
    expect([band.lo, band.hi].map((v) => v.toFixed(3)).join(','), '0.333,3.000');
    const cv50 = patientsRange(50);
    expect([cv50.lo, cv50.hi].map((v) => v.toFixed(2)).join(','), '0.41,1.95');
    expect([exposureRange(1.94, 50), exposureRange(0.42, 50)].join(','), ',');
    expect([exposureRange(1.96, 50)?.dir, exposureRange(0.40, 50)?.dir, exposureRange(3.1, null)?.dir].join(','),
      'above,below,above');
  });

  test('a combination names its drugs\' sources; a weekly multiple says so', async () => {
    const combo = study('CT26-X', [['Cpd-A 10 mg/kg QD', 600], ['Cpd-A 10 mg/kg QD + Cpd-B 5 mg/kg', 300]],
      [patients('Cpd-A')]);
    const r = compareRows(combo, 'CT26').find((g) => g.group.includes('+'))!;
    expect(r.source, 'Cpd-A: Projected');
    expect(compareFooter(compareRows(combo, 'CT26')).split('Projected').length, 2, 'the source said once');
    expect(shownMultiple({ratio: 5.3, tip: '', range: 'above', weekly: 'Mouse daily AUC…'}),
      '5.3× ↑, averaged over the week');
  });

  test('a combination\'s multiple: the same in the report as on screen, from its own group\'s AUC', async () => {
    const [one, both] = ['Cpd-A 10 mg/kg QD', 'Cpd-A 10 mg/kg QD + Cpd-B 5 mg/kg'];
    const combo = study('CT26-X', [[one, 600], [both, 300]], [patients('Cpd-A')], {table: 'pk', pk: new Map(),
      arms: new Map([[one, {exposure: x(10, 10)}]]),
      drugs: new Map([[both, [{key: drugKey('Cpd-A'), drug: 'Cpd-A', x: {exposure: x(10, 8)}, pk: null}]]])});
    const rows = compareRows(combo, 'CT26');
    const r = rows.find((g) => g.group === both)!;
    const onScreen = r.drugs!.find((d) => d.drug === 'Cpd-A')!.m!;
    expect(onScreen.ratio !== undefined, true);
    const s = compareSection(rows, ['CT26'])!;
    const line = s.rows.find((l) => l.some((c) => c.endsWith('(Cpd-A)') && c.includes('+')))!;
    expect(line !== undefined, true, s.rows.map((l) => l.join(' | ')).join(' / '));
    // (the figure itself: the report spells the arrow out)
    const figure = shownMultiple({ratio: onScreen.ratio!, tip: ''});
    expect(line[line.length - 2].startsWith(figure), true, line.join(' | '));
    expect(line[line.length - 1], 'Projected', 'its drug\'s source, not "Not given"');
  });

  test('two studies of one name: told apart by their IDs; a running one says so', async () => {
    expect(studyLabels([{name: 'CT26', id: 'PMX-1', stage: 2, day: 9}, {name: 'CT26', id: 'PMX-2', stage: 3,
      day: 16}, {name: 'MC38', id: 'PMX-3', stage: 3, day: 18}]).join('|'),
    'CT26 (PMX-1, running, day 9)|CT26 (PMX-2)|MC38');
    expect(studyLabels([{name: 'A', id: '', stage: 3, day: 1}, {name: 'A', id: '', stage: 3, day: 1}]).join('|'),
      'A|A #2');
  });

  test('a warning when the studies are of different models', async () => {
    expect(modelWarning(['HCT116', 'HCT116']), null);
    expect(modelWarning(['HCT116', 'B16F10'])!.startsWith('Different models'), true);
  });

  test('one hint per drug without the patients\' exposure', async () => {
    expect(JSON.stringify(hintRows(rows)), JSON.stringify([{drug: 'Cpd-B', key: 'cpdb', study: 'CRC-021'}]));
  });

  test('the sentence: from the numbers, by drug, never a judgement', async () => {
    const s = compareSentence(rows)!;
    expect(s.startsWith('Cpd-A gave a TGI of 40% at 10 mg/kg ('), true, s);
    expect(s.includes('60% at 30 mg/kg'), true, s);
    expect(s.endsWith('Paclitaxel has no mouse exposure; Cpd-B has no patients\' exposure yet.'), true, s);
    // a drug in two studies: each group's study named
    const again = study('CRC-030', [['Cpd-A 30 mg/kg QD', 500]], [patients('Cpd-A')], {table: 'pk', pk: new Map(),
      arms: new Map([['Cpd-A 30 mg/kg QD', {exposure: x(30, 40)}]])});
    const two = compareSentence([...rows, ...compareRows(again, 'HCT116')])!;
    expect(two.includes('at 10 mg/kg in CRC-017') && two.includes('in CRC-030'), true, two);
    const footer = compareFooter(rows);
    expect(footer.startsWith('Exposures are total, not free'), true, footer);
    for (const t of [s, footer])
      expect(/target|efficacious|clinical dose|margin|predict|did more|rose|fell/i.test(t), false, t);

    expect(compareSentence([]), null);
  });

  test('from the Overview: the studies read from their files', async () => {
    // a study file as Monitor saves it, and its entry as the library lists it: no text
    const file = (name: string, group: string) => JSON.stringify({format: 'pmx-study', version: 1, name, plan: null,
      looks: '', humane: null, log: [], meta: {...emptyMeta(), id: name}, rows: ([['Vehicle', 1000], [group, 500]] as
        [string, number][]).flatMap(([g, v]) => [1, 2, 3].map((i) => ({animal: `${g}-${i}`, group: g, day: 14,
        volume: v})))});
    const files = new Map([['a.pmxstudy', file('CRC-030', 'Cpd-A 10 mg/kg QD')],
      ['b.pmxstudy', file('CRC-031', 'Cpd-B 25 mg/kg')]]);
    const entries = [...files].map(([path, text]) => ({...entryOf(text, 1, {fileName: path} as DG.FileInfo)!,
      text: ''}));
    const read = async (f: DG.FileInfo) => files.get(f.fileName)!;
    expect([...(await compareTexts(entries, read)).texts.values()].join() === [...files.values()].join(), true);
    const v = await compareView(entries, read);
    try {
      const text = v.root.textContent ?? '';
      expect(text.includes('CRC-030') && text.includes('CRC-031') && text.includes('50% *'), true, text);
    } finally {
      v.close();
    }
    // a file deleted or damaged: that study left out, by name, the others shown
    const gone = async (f: DG.FileInfo) => f.fileName === 'b.pmxstudy' ? Promise.reject(new Error('gone')) : read(f);
    const some = await compareTexts(entries, gone);
    expect([some.texts.size, some.skipped.join()].join('|'), '1|CRC-031');
    const w = await compareView(entries, gone);
    try {
      expect((w.root.textContent ?? '').includes('CRC-031 could not be read'), true);
    } finally {
      w.close();
    }
  });

  test('the patients\' exposure added in Compare reaches the study open now', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['A1'])]);
    setMetaOn(df, {...emptyMeta(), id: 'CRC-021'});
    df.setTag(EXPOSURE_TAG, JSON.stringify({nca: {table: 'pk', arms: {}, pk: {}}}));
    const other = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['A1'])]);
    df.name = 'CRC-021';
    // a copy of the table: same ID, another name; left alone
    const copy = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['A1'])]);
    setMetaOn(copy, {...emptyMeta(), id: 'CRC-021'});
    copy.setTag(EXPOSURE_TAG, JSON.stringify({nca: {table: 'pk', arms: {}, pk: {}}}));
    copy.name = 'CRC-021 (2)';
    humansToOpen('CRC-021', 'CRC-021', new Map([patients('Cpd-B')]), [df, other, copy]);
    const tag = JSON.parse(df.getTag(EXPOSURE_TAG)!);
    expect([tag.nca.table, tag.humans.cpdb.auc].join(), 'pk,20000');
    expect(other.getTag(EXPOSURE_TAG) ?? null, null);
    expect(JSON.parse(copy.getTag(EXPOSURE_TAG)!).humans === undefined, true);
  });

  test('the patients\' exposure saved on the file as it is now, not as Compare read it', async () => {
    // the file saved since: a new NCA table and a new note; only the patients' exposure changes
    const now = JSON.stringify({format: 'pmx-study', name: 'CRC-021', notes: 'saved since', exposure:
      {nca: {table: 'pk-new', arms: {}, pk: {}}}});
    const out = JSON.parse(withHumans(now, new Map([patients('Cpd-B')])));
    expect([out.notes, out.exposure.nca.table, out.exposure.humans.cpdb.auc].join('|'), 'saved since|pk-new|20000');
    expect(JSON.parse(withHumans(now, new Map())).exposure.humans === undefined, true);
  });
});

category('Compare: growth, doses and controls against the atlas', () => {
  test('a control against the atlas: within its range, or how many times faster or slower', async () => {
    // the atlas: log volume growing 0.15 to 0.30 a day (doubling 2.3 to 4.6 d), median 0.2 (3.5 d)
    const atlas = {q025: 0.15, median: 0.2, q975: 0.3};
    const within = controlVsAtlas(0.2, atlas);
    expect([within.side, sideText(within), rangeText(within)].join('|'), 'within|Within the atlas range|2.3–4.6 d');
    const fast = controlVsAtlas(0.4, atlas);
    expect([fast.side, sideText(fast)].join('|'), 'faster|2× faster');
    const slow = controlVsAtlas(0.1, atlas);
    expect([slow.side, sideText(slow)].join('|'), 'slower|2× slower');
    // no prediction (no model, the atlas failed) or a control that did not grow: nothing said against the atlas
    expect(controlVsAtlas(0.2, null).side, null);
    // a control that did not grow: the furthest below the range, said so
    const flat = controlVsAtlas(-0.01, atlas);
    expect([flat.side, flat.doubling, sideText(flat)].join('|'), 'slower||Did not grow');
    expect(driftSentence([{study: 'C', model: 'x', check: flat}])!.startsWith('C: its control did not grow'), true);
    // one plain sentence, only for the controls that drifted
    expect(driftSentence([{study: 'A', model: 'x', check: within}]), null);
    expect(driftSentence([{study: 'B', model: 'x', check: fast}])!.startsWith('B: its control doubled every 1.7 d, ' +
      'faster than the atlas predicts'), true);
  });

  test('the control’s growth: the mean of its animals’ rates, its last day and its day-0 median', async () => {
    const s: PreviousStudy = {name: 'S', line: null, meta: null, rows: [0, 5, 10].flatMap((d) => [
      {animal: 'c1', group: 'Vehicle', day: d, volume: 100 * Math.exp(0.2 * d)},
      {animal: 'c2', group: 'Vehicle', day: d, volume: 120 * Math.exp(0.1 * d)},
      {animal: 't1', group: 'Cpd-A 10 mg/kg', day: d, volume: 110}])};
    const g = controlGrowth(s);
    expect([sig(g.slope!), g.lastDay, g.vRand, g.n].join('|'), '0.15|10|110|2');
    // read to the evaluation day, as Analyze reads it: before day 0 and past the endpoint (where a control slows)
    // left out
    const slowing: PreviousStudy = {name: 'S', line: null, meta: null, rows: [[-3, 60], [0, 100], [7, 400], [14, 1600],
      [21, 1700]].map(([day, volume]) => ({animal: 'c1', group: 'Vehicle', day, volume}))};
    const cut = controlGrowth(slowing);
    expect([cut.toDay, sig(cut.slope!)].join('|'), `14|${sig(Math.log(16) / 14)}`);
    expect(controlGrowth(slowing, 1000, 7).toDay, 7);
  });

  test('a control read as Analyze read it, when the study kept its analysis', async () => {
    const atlas = {q025: 0.15, median: 0.2, q975: 0.3};
    // its rows alone (a plateau after the endpoint) would read slower than the atlas
    const rows = [[0, 100], [7, 140], [14, 150]].map(([day, volume]) => ({animal: 'c1', group: 'Vehicle', day,
      volume}));
    const plain: PreviousStudy = {name: 'S', line: null, meta: null, rows};
    expect(controlCheck(plain, atlas).check.side, 'slower');
    // kept: the analysis's rate and range, classified as Analyze's verdict
    const kept = {slope: 0.18, q025: 0.12, median: 0.17, q975: 0.25, endpoint: 1000};
    const s: PreviousStudy = {...plain, exposure: {effect: {day: 14, tgi: {}, control: kept}}};
    const {check, fromAnalysis} = controlCheck(s, atlas);
    expect([fromAnalysis, check.side].join('|'), 'true|within');
    expect(slopeSide(kept.slope, kept), 'As predicted');
    for (const x of [0.05, 0.4]) {
      const k = {...s, exposure: {effect: {day: 14, tgi: {}, control: {...kept, slope: x}}}};
      expect(controlCheck(k, atlas).check.side, {'As predicted': 'within', 'Faster than predicted': 'faster',
        'Slower than predicted': 'slower'}[slopeSide(x, kept)]);
    }
  });

  test('doses: each drug alone by its unit, the lowest first; combinations and groups without a dose left out',
    async () => {
      const rows = [...compareRows(study('S1', [['Cpd-A 30 mg/kg QD', 400], ['Cpd-A 10 mg/kg QD', 600],
        ['Cpd-A + Cpd-B', 300]]), 'x'), ...compareRows(study('S2', [['Cpd-A 3 mg/kg QD', 900]]), 'x')];
      const groups = doseGroups(rows);
      expect(groups.map((g) => g.name).join(), 'Cpd-A (mg/kg)');
      expect(groups[0].rows.map((r) => `${r.doseValue}:${r.study}`).join(), '3:S2,10:S1,30:S1');
    });

  test('growth: a series per group and study, the treatment’s colour, the second study dashed; V/V0',
    async () => {
      const one = study('S1', [['Cpd-A 10 mg/kg QD', 500]]);
      const two = study('S2', [['Cpd-A 10 mg/kg QD', 700]]);
      const g = growthFrame([{s: one, label: 'S1'}, {s: two, label: 'S2'}], false)!;
      const series = [...new Set(g.df.col('series')!.toList())];
      expect(series.length, 4);
      expect(g.dashed.every((x) => x.endsWith('· S2')) && g.dashed.length === 2, true, g.dashed.join());
      const a = series.filter((x) => !x.startsWith('Control'));
      expect(g.colors[a[0]] === g.colors[a[1]], true, 'one treatment, one colour');
      // a group the study recoloured: drawn in its study's colour, as Analyze draws it
      const own = {...one, meta: {...emptyMeta(), colors: {'Cpd-A 10 mg/kg QD': '#ff7f0e'}}};
      const c = growthFrame([{s: own, label: 'S1'}], false)!.colors;
      expect(c[[...Object.keys(c)].find((x) => x.startsWith('Cpd-A'))!], '#ff7f0e', JSON.stringify(c));
      // relative: each animal over its first volume (one day each here: 1)
      expect(growthFrame([{s: one, label: 'S1'}], true)!.df.col('median')!.toList().every((v) => v === 1), true);
    });
});

category('Report: exposure and patients', () => {
  const ok = (s: ExposureSection) => [s.sentence ?? '', ...s.rows.flat(), ...s.sources, ...s.notes.slice(1)]
    .every((t) => !/target|efficacious|clinical dose|margin|predict|did more|rose|fell/i.test(t));

  test('no exposure: no section', async () => {
    expect(exposureSection({sentence: null, curve: false, rows: [{group: 'Cpd-A 10 mg/kg QD', tgi: 40,
      x: undefined, multiple: null, human: undefined}]}), null);
  });

  test('Analyze: the groups, a projected value, the curve, the sources', async () => {
    const s = exposureSection({sentence: 'At an exposure within the patients\' range (1.4×), 30 mg/kg gave a TGI ' +
      'of 62 %.', curve: true, table: 'pk_conc', rows: [{group: 'Cpd-A 30 mg/kg QD', tgi: 62,
      x: {exposure: x(30, 40)}, multiple: {ratio: 4.6, tip: '', range: 'above'},
      human: patients('Cpd-A')[1]}]})!;
    expect(s.head.join('|'), 'Group|TGI from the analysis|Mouse AUC per dose|× patients|' +
      'Patients\' exposure from');
    expect(s.rows[0].join('|'), 'Cpd-A 30 mg/kg QD|62%|40 mg/L·h (AUCinf)|4.6× ↑ above the patients\' range|' +
      'Projected');
    expect(s.curve, 'Within the tested doses only.');
    expect(s.notes[0].startsWith('The mice and the patients are set side by side'), true);
    expect(s.notes.some((n) => n.startsWith('Exposures are total, not free')), true);
    expect(s.notes.some((n) => n.includes('Cpd-A is a projection')), true);
    expect(s.notes.some((n) => n.includes('between the doses tested')), true);
    expect(s.sources.join(' '), 'Mouse exposure: NCA, table pk_conc. Patients\' exposure of Cpd-A: Projected, ' +
      'CV 45%.');
    expect(ok(s), true, JSON.stringify(s));
  });

  test('a per-cycle reference row, and a group not placed: its reason', async () => {
    const row: ReferenceRow = {drug: 'trastuzumab', synonyms: ['Herceptin'], regimen: '6 mg/kg IV every 3 weeks',
      auc: 43056000,
      aucWindow: 'AUCss,0-21 days', cycleDays: 21, aucCv: null, cmax: null, cmaxCv: null, fu: null, matrix: 'plasma',
      stat: 'median', tier: 'label', source: 'FDA label 12.3'};
    const h = {...fromReference('Herceptin', row), reference: {row, edited: false}};
    const s = exposureSection({sentence: null, curve: false, table: 'pk', rows: [{group: 'Herceptin 10 mg/kg QW',
      tgi: 50, x: {exposure: {...x(10, 7000), drug: 'Herceptin'}}, multiple: {ratio: 1.1, tip: ''}, human: h},
    {group: 'Herceptin 30 mg/kg', tgi: 70, x: undefined, multiple: {why: 'its schedule is not read'}, human: h}]})!;
    expect(s.rows[0][4], 'FDA label 12.3, 6 mg/kg IV every 3 weeks');
    expect(s.rows[1][3], '— (its schedule is not read)');
    expect(s.notes.some((n) => n.includes('projection')), false);
    expect(s.curve, null);
    expect(ok(s), true, JSON.stringify(s));
  });

  test('the patients\' exposure without NCA: no mouse exposure named as a source', async () => {
    const s = exposureSection({sentence: null, curve: false, table: null, rows: [{group: 'Cpd-A 30 mg/kg QD', tgi: 62,
      x: undefined, multiple: {why: 'no NCA'}, human: patients('Cpd-A')[1]}]})!;
    expect(s.sources.join(' '), 'Patients\' exposure of Cpd-A: Projected, CV 45%.');
  });

  test('Compare: the studies, the models\' warning', async () => {
    const a = study('CRC-017', [['Cpd-A 30 mg/kg QD', 400]], [patients('Cpd-A')], {table: 'pk', pk: new Map(),
      arms: new Map([['Cpd-A 30 mg/kg QD', {exposure: x(30, 40)}]])});
    const s = compareSection(compareRows(a, 'HCT116'), ['HCT116', 'B16F10'])!;
    expect(s.head.slice(0, 3).join('|'), 'Study|Model|Group');
    expect(s.rows[0].slice(0, 2).join('|'), 'CRC-017|HCT116');
    // the medians' TGI, starred and said so; each study's value named with its study
    expect(s.rows[0][3], '60% *');
    expect(s.notes.includes('* The study kept no analysis: 1 − T/C of the median volumes on day 14.'), true,
      s.notes.join(' | '));
    expect(s.sources[1], 'Patients\' exposure of Cpd-A (CRC-017): Projected, CV 45%.');
    expect(s.warning!.startsWith('Different models'), true);
    expect(s.sources[0], 'Mouse exposure: NCA kept with each study.');
    expect(ok(s), true, JSON.stringify(s));
  });
});

category('Analyze: a group\'s translation details', () => {
  const base = {group: 'Cpd-A 30 mg/kg QD', read: 'Cpd-A · 30 mg/kg · QD', tgi: 62, day: 21, x: undefined,
    table: 'pk_conc', live: true, multiple: null, human: undefined, curve: null};
  const headings = (d: GroupDetails) => d.blocks.map((b) => b.heading).join('|');

  test('the answer’s p-value line says the side its result says; a weight gained is signed', async () => {
    const sub = (result: string, p: number) => groupDetails({...base, result, p}).verdict?.sub ?? '';
    expect(sub('Slows growth', 0.001).startsWith('Clearly below the control'), true);
    expect(sub('Grows faster', 0.001).startsWith('Clearly above the control'), true);
    expect(sub('No clear effect', 0.4).startsWith('Not clearly different from the control'), true);
    const weight = (nadir: number) => groupDetails({...base, nadir}).facts.find((f) => f.label === 'Weight')?.value;
    expect(weight(2), '+2% at its lowest');
    // its treatment in plain words under its name
    expect(plainArm('Cpd-A · 30 mg/kg · QD · p.o.'), 'Cpd-A · 30 mg/kg · once a day · by mouth');
    expect(plainArm('Cpd-A · 5 mg/kg · Q3D · i.p. + anti-PD-1 · 10 mg/kg · BIW'),
      'Cpd-A · 5 mg/kg · every 3 days · into the abdomen + anti-PD-1 · 10 mg/kg · twice a week');
    expect(weight(-6), '−6% at its lowest');
  });

  test('full: the facts on top, both exposures, the multiple, the curve folded', async () => {
    const d = groupDetails({...base, x: {exposure: x(30, 40)}, human: patients('Cpd-A')[1],
      multiple: {ratio: 2, tip: 'Mouse daily AUC: AUCinf per dose × 1 dose a day (QD).', range: 'above'},
      curve: 'between'});
    expect(headings(d), 'Mouse exposure|Patients\' exposure|Against the patients|Descriptive curve');
    expect(d.facts.map((f) => `${f.label}: ${f.value}`).join('|'), 'TGI: 62% on day 21');
    const lines = d.blocks.flatMap((b) => b.lines);
    expect(lines.includes('AUCinf: 40 mg/L·h') && lines.includes('AUClast: 36 mg/L·h'), true, lines.join(' | '));
    expect(lines.includes('Mouse daily AUC: AUCinf per dose × 1 dose a day (QD).'), true, 'which AUC: the multiple\'s');
    expect(lines.includes('NCA, table pk_conc'), true);
    expect(lines.includes('Source: Projected') && lines.includes('CV 45%'), true);
    expect(lines.includes('2× ↑ above the patients\' range'), true);
    expect(lines.every((l) => !/target|efficacious|clinical dose|margin|predict|did more|rose|fell/i.test(l)), true);
  });

  test('no NCA: what the group has, the reason (what would add to it is said once, under the groups)', async () => {
    const d = groupDetails({...base, x: {why: 'no group at 30 mg/kg'}, human: patients('Cpd-A')[1]});
    expect(headings(d), 'Mouse exposure|Patients\' exposure');
    expect(d.blocks[0].lines[0], 'no group at 30 mg/kg');
  });

  test('no patients\' exposure', async () => {
    const d = groupDetails({...base, x: {exposure: {...x(30, 40), aucInf: null}}, live: false});
    expect(headings(d), 'Mouse exposure');
    expect(d.blocks[0].lines.some((l) => l.startsWith('AUCinf')), false);
    expect(d.blocks[0].lines.includes('NCA, saved with the study (table pk_conc)'), true);
  });

  test('a combination: its reading and why it has no multiple', async () => {
    const d = groupDetails({...base, group: 'Cpd-A 30 mg/kg + Cpd-B 10 mg/kg',
      read: 'Cpd-A · 30 mg/kg + Cpd-B · 10 mg/kg',
      multiple: {why: 'a combination: its drugs\' exposures are not added up into one multiple'}, curve: 'off'});
    expect(headings(d), 'Against the patients|Descriptive curve');
    expect(d.read.includes(' + '), true);
    expect(d.blocks[1].lines[0].startsWith('Not on the curve'), true);
  });

  test('the panel: name and colour, its answer first, plain facts, data only (no formula, no link)', async () => {
    const d = groupDetails({...base, group: 'LD4172 20 mg/kg i.p. QD + anti-PD-1',
      read: 'LD4172 · 20 mg/kg + anti-PD-1',
      color: '#1f77b4', n: 8, tgi: 108, day: 6, tc: {value: -8, lower: -20, upper: 4}, result: 'Slows growth',
      step: {words: 'More than its two drugs acting independently', why: 'Bliss: why'}, responders: '5 of 8',
      x: {exposure: x(20, 40)}, tgiCi: {lower: 95, upper: 121}, nadir: -12.4, p: 0.0004,
      removed: '2 removed on day 7 (weight)'});
    expect(d.facts.map((f) => f.label).join('|'), 'TGI|Size vs control|Against its drugs|Responders|Animals|Weight');
    // the combination step's own why on its fact's hover only
    expect(d.facts.find((f) => f.label === 'Against its drugs')?.tip, 'Bliss: why');
    const el = detailsElement(d, [{title: 'Animals', content: () => document.createElement('div')}]);
    const title = el.querySelector('.pmx-object-title');
    expect(title?.textContent, 'LD4172 20 mg/kg i.p. QD + anti-PD-1');
    expect(title?.querySelector('span[style*="background"]') !== null, true, 'the colour dot');
    const facts = el.querySelector('.pmx-object-facts')?.textContent ?? '';
    for (const s of ['Animals8 · 2 removed on day 7 (weight)', 'TGI108% (95–121%) on day 6', 'Weight−12% at its lowest',
      '5 of 8 measured on day 6'])
      expect(facts.includes(s), true, `${s} in ${facts}`);
    // the answer before the facts, in the groups table's words, with its p-value
    const verdict = el.querySelector('.pmx-object-verdict');
    expect(verdict?.firstChild?.textContent, 'Slows growth: TGI 108% on day 6');
    expect(verdict?.textContent?.includes('Clearly below the control (p < 0.001)'), true, verdict?.textContent ?? '');
    expect(!!(verdict!.compareDocumentPosition(el.querySelector('.pmx-object-facts')!) &
      Node.DOCUMENT_POSITION_FOLLOWING), true, 'the answer first');
    // no formula and no link anywhere: the panel is data
    expect(/ΔT|ΔC|Method|exposure…|computed/.test(el.textContent ?? ''), false, el.textContent ?? '');
    const panes = [...el.querySelectorAll('.d4-accordion-pane-header')].map((h) => h.textContent?.trim());
    expect(panes.join('|'), 'Animals|Exposure');
    // nothing to say: no sections at all
    const bare = groupDetails({...base, tgi: null});
    expect(detailsElement(bare).querySelector('.pmx-object-sections') === null, true);
  });
});

category('Context panel: an animal', () => {
  const study = () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M1', 'M2']),
      DG.Column.fromStrings('~group', ['Cpd-A', 'Cpd-A', 'Cpd-A', 'Vehicle']),
      DG.Column.fromList(DG.TYPE.INT, 'day', [7, 0, 3, 0]),
      DG.Column.fromList(DG.TYPE.FLOAT, 'volume', [400, 100, 900, 120]),
      DG.Column.fromList(DG.TYPE.BOOL, 'Excluded', [false, false, true, false]),
      DG.Column.fromStrings('Mark reason', ['', '', 'caliper slip', ''])]);
    return {study: df, animal: 'animal', group: '~group', day: 'day', volume: 'volume', bodyWeight: null};
  };

  test('its readings by day, the excluded one left out with its reason', async () => {
    const r = animalReadings(study(), 'M1');
    expect(r.days.join(','), '0,7');
    expect(r.volumes.join(','), '100,400');
    expect(r.excluded, 1);
    expect(r.reasons.join(), 'caliper slip');
    expect(r.group, 'Cpd-A');
  });

  test('blinded: no group and no colour', async () => {
    const el = animalPanel(study(), 'M1', {showGroup: false, color: '#1f77b4'});
    const text = el.textContent ?? '';
    expect(text.includes('Cpd-A'), false, text);
    expect(el.querySelector('span[style*="background"]') === null, true, 'no colour dot');
    expect(text.includes('400 mm³ on day 7'), true, text);
    // Monitor: no answer line, never made up
    expect(el.querySelector('.pmx-object-verdict') === null, true);
    expect(animalPanel(study(), 'M1', {showGroup: true}).textContent?.includes('Cpd-A'), true);
  });
});
