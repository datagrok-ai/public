import * as DG from 'datagrok-api/dg';
import {PmxDay} from '../pmx-objects';
import * as grok from 'datagrok-api/grok';
import {category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_RULES, findOutliers, rulesProblem, StudyColumns} from '../outlier-rules';
import {AnimalSelection, animalMarks, appendPasted, caliperVolume, emptyStudy, exclusionsOf, mark, mirrorSelection,
  dayFormat, findingsOf, pastedNumber, studyDataLine, withoutExcluded} from '../study-table';
import {EXCLUDED, FINISHED_TAG, FLAGGED, LINE_TAG, MARK_REASON, REVEALED, STUDY_TAG, WEIGHT_TAG, addLog, groupsShown,
  lineOf, logOf, freeName, parseStudy, roleColumns, studyOf, tableOf, analysisFieldsOf, setAnalysisOn}
  from '../study-file';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {exampleStudy} from './example-study';
import {filePreview, studySummary} from '../file-preview';
import {EXAMPLE, isBodyWeight} from '../ui';
import {_package} from '../package';
import {armDoses, armParts, dateOf, dosingDays, emptyMeta, metaOf, parseLine, parseMeta, setMetaOn, startFrom,
  weekdayOf, fromName, withGroups, withException, withGiven, doseForTotal, renamedDetails, StudyMeta, treatmentOf, compoundsOf,
  onlyGroups}
  from '../study-meta';
import {parseArmName} from '../arm-name';
import {dosedDays, dosesPerDay, excepted, groupWeights, readSchedule} from '../schedule';
import {entryOf, libraryOf, matches} from '../library';
import {animalOf, spiderLines} from '../spider-chart';
import {dayActions, dayToday, icsOf, startMismatch, weeksGrid} from '../calendar';
import {Scope} from '@datagrok-libraries/u2';
import {sheetHtml, sheetRows} from '../sheet';

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
  test('The example study opens: 160 rows, 4 groups, days 0 to 12, a body weight', async () => {
    const df = DG.DataFrame.fromCsv(await _package.files.readAsText(EXAMPLE.file));
    const {group, day, volume} = roleColumns(df);
    expect(df.rowCount, 160, 'rows');
    expect(group.categories.length, 4, 'groups');
    expect([...new Set(day.toList())].sort((a, b) => a - b).join(), '0,3,6,9,12', 'days from randomization');
    expect(volume.stats.missingValueCount, 0, 'no volume missing');
    expect(df.columns.names().filter(isBodyWeight).length, 1, 'one body weight column');
  });

  test('a study file previewed: its animals, groups and days; no group named while blinded; a newer file said so',
    async () => {
      const df = DG.DataFrame.fromCsv(await _package.files.readAsText(EXAMPLE.file));
      df.setTag(WEIGHT_TAG, 'body_weight');
      const saved = parseStudy(JSON.stringify(studyOf(df, 1500)));
      const lines = studySummary(saved);
      const all = lines.join(' | ');
      expect(lines.includes('32 animals in 4 groups, measured on 5 days, the last day 12'), true, all);
      expect(lines.includes('Vehicle: 8 animals'), true, all);
      expect(lines.includes('Body weight measured') && lines.includes('No plan attached'), true, all);
      // blinded: no group's name in any line, nor how many there are
      const names = [...new Set(saved.rows.map((r) => r.group!))];
      const blind = studySummary({...saved, meta: {...emptyMeta(), blinded: true}}).join(' | ');
      expect(names.some((g) => blind.includes(g)) || /\d+ groups?\b/.test(blind), false, blind);
      const revealed = studySummary({...saved, meta: {...emptyMeta(), blinded: true},
        log: [{when: '2026-11-12', who: 'x', what: REVEALED}]});
      expect(revealed.includes('Vehicle: 8 animals'), true, 'revealed, still unnamed');
      expect(studySummary({...saved, rows: []})[0], '0 animals, measured on 0 days');
      // blinded with a plan attached: the plan is said to be there, not named
      const planned = studySummary({...saved, plan: {...saved.plan, name: 'Secret plan'} as any,
        meta: {...emptyMeta(), blinded: true}});
      expect(planned.includes('A plan attached') && !planned.join(' ').includes('Secret plan'), true,
        planned.join(' | '));
      // a study that reads: its lines and the button, which opens it
      let opened = '';
      const ok = await filePreview(DG.FileInfo.fromString('S.pmxstudy', JSON.stringify(studyOf(df, 1500))), 'study',
        (t) => opened = t);
      expect(ok.root.textContent!.includes('Vehicle: 8 animals'), true, ok.root.textContent!);
      ok.root.querySelector('button')!.click();
      await new Promise((r) => setTimeout(r, 0));
      expect(opened.length > 0, true, 'the button opens the study');
      // the preview itself: a file that is no study, or damaged, says so in a sentence
      for (const [text, says] of [['{"format": "other"}', 'Not a PMX Studio study'], ['{', 'cannot be read']]) {
        const v = await filePreview(DG.FileInfo.fromString('x.pmxstudy', text), 'study', () => {});
        expect(v.root.textContent!.includes(says), true, v.root.textContent!);
      }
      let said = '';
      try {
        parseStudy(JSON.stringify({...saved, version: 99}));
      } catch (e: any) {
        said = e.message;
      }
      expect(said.includes('update the PMX Studio package'), true, said);
    });

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
    // values that look like headings are values: nothing dropped
    expect(appendPasted(df, 'Mouse 4\tGroup 2\t3\t160'), 1);
    expect(df.get('animal', 3), 'Mouse 4');
  });

  test('a pasted number is read or refused, never guessed: nothing is added when one cell cannot be read', async () => {
    expect(pastedNumber('101.5'), 101.5);
    expect(pastedNumber('101,5'), 101.5);
    expect(pastedNumber('1,234.5'), 1234.5);
    expect(pastedNumber('1.234,5'), 1234.5);
    expect(pastedNumber(' n/a '), null);
    expect(Number.isNaN(pastedNumber('1,234')), true, '1.234 or 1234: not guessed');
    expect(Number.isNaN(pastedNumber('120 mm3')), true);
    const df = emptyStudy();
    const refused = (text: string) => {
      try {
        appendPasted(df, text);
        return '';
      } catch (e: any) {
        return String(e.message);
      }
    };
    expect(refused('M1\tVehicle\t0\t100\nM2\tVehicle\t0\t1,234').includes('row 2, volume: "1,234"'), true);
    expect(refused('animal\t0\t3\t7\nM1\t100\t150\t220').includes('wide table'), true);
    expect(refused('M1;Vehicle;0;100').includes('no tabs'), true);
    expect(df.rowCount, 0, 'none of the three added a row');
  });

  test('the jump rule reads a change against the group’s over the same days; excluded rows are out of it', async () => {
    // measured weekly, every tumor about four times larger: growth, not a typo; M3 tenfold is one
    const s: StudyColumns = {animal: [], group: [], day: [], volume: [], weight: null};
    for (const [a, v7] of [['M1', 400], ['M2', 420], ['M3', 1000], ['M4', 380]] as [string, number][]) {
      for (const [day, v] of [[0, 100], [7, v7]]) {
        s.animal.push(a);
        s.group.push('Vehicle');
        s.day.push(day);
        s.volume.push(v);
      }
    }
    const jumps = (x: StudyColumns) => findOutliers(x, {...DEFAULT_RULES, group: {...DEFAULT_RULES.group, on: false}});
    expect(jumps(s).length, 0, 'a fourfold week is what the group did; tenfold is under three times that');
    s.volume[5] = 4000; // M3's day 7
    expect(jumps(s).map((f) => f.row).join(','), '5');
    // an animal alone (no group to compare with) is read against no change, as before
    expect(jumps({animal: ['M1', 'M1'], group: ['A', 'A'], day: [0, 3], volume: [100, 320], weight: null}).length, 1);
    // an animal out of the analysis in part is said to be so, not "Excluded"
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M2', 'M2', 'M3']),
      DG.Column.fromList('bool', EXCLUDED, [true, true, true, false, false]),
      DG.Column.fromList('bool', FLAGGED, [false, false, false, false, true])]);
    const marks = animalMarks({study: df, animal: 'animal', group: 'animal', day: 'animal', volume: 'animal',
      bodyWeight: null});
    expect([marks.get('M1'), marks.get('M2'), marks.get('M3')].join(','), 'excluded,partly,flagged');
  });

  test('a study with its own column names and extra columns saves by role and keeps them', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('Mouse', ['M1', 'M2']),
      DG.Column.fromStrings('Arm', ['A', 'B']), DG.Column.fromList('double', 'Day', [0, 0]),
      DG.Column.fromList('double', 'TV', [100, 110]), DG.Column.fromList('double', 'length', [6.5, 6.25]),
      DG.Column.fromStrings('cage', ['C1', 'C2'])]);
    df.name = 'S';
    let refused = false;
    try {
      studyOf(df, 2000);
    } catch {
      refused = true;
    }
    expect(refused, true, 'no roles yet: the usual names are not there');
    df.setTag(STUDY_TAG, JSON.stringify({animal: 'Mouse', group: 'Arm', day: 'Day', volume: 'TV', bodyWeight: null}));
    const back = tableOf(parseStudy(JSON.stringify(studyOf(df, 2000))));
    expect(back.get('animal', 1), 'M2');
    expect(back.get('volume', 1), 110);
    expect(back.get('length', 0), 6.5, 'the caliper reading travels with the study');
    expect(back.get('cage', 1), 'C2');
    let damaged = false;
    try {
      parseStudy(JSON.stringify({...studyOf(df, 2000), log: [{who: 'x'}]}));
    } catch {
      damaged = true;
    }
    expect(damaged, true, 'a log line with no time is refused on opening, not deep in a view');
  });

  test('excluded rows leave the engine’s inputs and are reported with their reasons', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M2', 'M2']),
      DG.Column.fromStrings('group', ['A', 'A', 'A', 'A']), DG.Column.fromList('double', 'day', [0, 3, 0, 3]),
      DG.Column.fromList('double', 'volume', [100, 150, 90, 900]),
      // a study file's own Excluded column, as a lab writes it
      DG.Column.fromStrings(EXCLUDED, ['x', 'x', '', 'yes']),
      DG.Column.fromStrings(MARK_REASON, ['dosing error', 'dosing error', '', 'misread'])]);
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

  test('marks are ticks in visible columns; flagging keeps an exclusion; clearing removes it; all logged', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M1'])]);
    mark(df, [1], 'excluded', 'misread', 'M1 on day 3');
    mark(df, [0, 1, 2], 'flagged', 'thin', 'animal M1');
    expect(df.col(EXCLUDED)!.type, DG.TYPE.BOOL);
    expect(df.col(EXCLUDED)!.toList().join(','), 'false,true,false');
    expect(df.col(FLAGGED)!.toList().join(','), 'true,true,true');
    expect(df.get(MARK_REASON, 1), 'misread', 'the exclusion keeps its reason');
    expect(df.get(MARK_REASON, 0), 'thin');
    expect(withoutExcluded({study: df}).study.rowCount, 2);
    mark(df, [0, 1, 2], null, '', 'animal M1');
    expect([...df.col(EXCLUDED)!.toList(), ...df.col(FLAGGED)!.toList()].every((v) => !v), true);
    expect(df.rowCount, 3, 'nothing deleted');
    expect(JSON.parse(df.getTag('pmx.log') || '[]').length, 3);
  });

  test('the spider plot: each animal from its volume at randomization, excluded rows left out', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1', 'M1', 'M2', 'M2']),
      DG.Column.fromStrings('group', ['A', 'A', 'A', 'B', 'B']), DG.Column.fromList('double', 'day', [7, 0, 3, 0, 3]),
      DG.Column.fromList('double', 'volume', [50, 100, 150, 0, 80]),
      DG.Column.fromList('bool', EXCLUDED, [false, false, true, false, false])]);
    const lines = spiderLines({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume',
      bodyWeight: null});
    expect(lines.length, 1, 'M2 starts at 0 mm³: no line');
    expect(lines[0].days.join(','), '0,7', 'by day, day 3 excluded');
    expect(lines[0].change.join(','), '0,-50');
    // a click on a point, or on the line between two (no point): the animal either way; a bar's needs its point
    expect(animalOf({seriesType: 'line', seriesId: 'M1', data: {animal: 'M1'}}), 'M1');
    expect(animalOf({seriesType: 'line', seriesId: 'M1'}), 'M1', 'the line between points, hovered');
    expect(animalOf({seriesType: 'line', seriesIndex: 1}, ['M1', 'M2']), 'M2', 'the line between points, clicked');
    expect(animalOf({seriesType: 'bar', seriesId: 'series0'}) === undefined, true, 'a bar without its point');
  });

  test('the randomization sheet: each animal once, at randomization, by group', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M2', 'M2', 'M10', 'M1', 'M3']),
      DG.Column.fromStrings('group', ['B', 'B', 'A', 'A', 'B']), DG.Column.fromList('double', 'day', [-2, 0, 3, 0, 4]),
      DG.Column.fromList('double', 'volume', [80, 101, 150, 99, 120])]);
    const rows = sheetRows(df);
    expect(rows.map((r) => `${r.group}:${r.animal}:${r.volume}`).join(','), 'A:M1:99,A:M10:150,B:M2:101,B:M3:120',
      'day 0 over day -2; an animal first seen later keeps its first');
    const html = sheetHtml('S <1>', {...emptyMeta(), id: 'PMX-2026-002', start: '2026-11-09', groups:
      [{group: 'A', compound: 'anti-PD-1', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'Q3D×4'}]}, df, '42');
    expect(html.includes('S &lt;1&gt;'), true, 'names escaped');
    expect(html.includes('anti-PD-1 · 10 mg/kg · i.p. · Q3D×4'), true);
    expect(html.includes('Seed 42'), true);
  });

  test('caliper readings make the volume, live', async () => {
    const df = emptyStudy();
    df.rows.addNew(['M1', 'A', 0, null, 20], false);
    expect(await caliperVolume(df, '(add the column)', '(add the column)'), 'volume', 'the empty volume replaced');
    df.set('length', 0, 10);
    df.set('width', 0, 6);
    await new Promise((r) => setTimeout(r, 50));
    expect(df.get('volume', 0), 180, '10 × 6² / 2');
    const kept = DG.DataFrame.fromColumns([DG.Column.fromList('double', 'volume', [5]),
      DG.Column.fromList('double', 'L', [4]), DG.Column.fromList('double', 'W', [2])]);
    expect(await caliperVolume(kept, 'L', 'W'), 'volume (caliper)', 'a volume with values stays');
    expect(kept.get('volume (caliper)', 0), 8);
  });

  test('the calendar’s dates, doses and notes', async () => {
    const g = (group: string, schedule: string) => ({group, compound: '', dose: null, unit: '', route: '', schedule});
    const doses = dosingDays([g('A', 'Q3D×4'), g('B', 'qw x 3'), g('C', 'BIW'), g('D', 'QD×2')], 14);
    expect([...doses.keys()].sort((a, b) => a - b).join(','), '0,1,3,6,7,9,10,14');
    expect(doses.get(0)!.join(','), 'A,B,C,D', 'every group on day 0, BIW too');
    expect(doses.get(9)!.join(','), 'A');
    expect(dosingDays([g('A', 'Q7D×10')], 20).size, 3, 'nothing past the last day');
    expect(dateOf('2026-11-09', 22)!.toISOString().slice(0, 10), '2026-12-01');
    expect(dateOf(null, 3), null);
    expect(startFrom('2026-11-12', 3), '2026-11-09');
    expect(weekdayOf('2026-11-09'), 0, 'a Monday');
    expect(dayToday('2026-11-09', new Date(2026, 10, 21, 23, 30)), 12, 'today, by the local date');
    expect(dayToday('2026-11-09', new Date(2026, 10, 2)), -7, 'before day 0');
    expect(dayToday(null), null);
    // QA pass 4: a finished table measured to day 21 by 25 Nov, day 0 dated 30 Nov: the date is not this table's
    const nov25 = new Date(2026, 10, 25);
    expect(startMismatch('2026-11-30', 21, nov25)?.includes('measurements up to day 21'), true);
    expect(startMismatch('2026-11-04', 21, nov25), null, 'day 21 is today');
    expect(startMismatch('2026-11-30', null, nov25), null, 'no measurements');
    expect(startMismatch(null, 21, nov25), null, 'no date');
    const ics = icsOf('S', {...emptyMeta(), start: '2026-11-09', groups: [g('A', 'QD×2')]}, 3, [[0, 'Randomize']]);
    expect(ics.split('BEGIN:VEVENT').length - 1, 3, 'randomize and two doses, all-day');
    expect(ics.includes('DTSTART;VALUE=DATE:20261110'), true);
    // past the expected end an event says it is conditional; its UID stays the same, so a later export updates it
    const past = icsOf('S', {...emptyMeta(), start: '2026-11-09', groups: [g('A', 'QD×2')]}, 3, [[0, 'Randomize']],
      0);
    expect(past.includes('SUMMARY:S: Dose: A (day 1)\\, only if the study runs past its expected end (day 0)'), true);
    expect(past.includes('SUMMARY:S: Dose: A (day 0)\r\n'), true, 'the expected end\'s own day is not conditional');
    expect(past.split('\r\n').filter((l) => l.startsWith('UID:')).join(),
      ics.split('\r\n').filter((l) => l.startsWith('UID:')).join(), 'the same UIDs');
    const notes = parseMeta({notes: [{day: 3, kind: 'Incident', animal: 'M4', text: 'ulcerated', who: 'x', when: ''},
      {day: 'x', text: 'no day'}, {day: 5, kind: 'Party', text: 'odd kind'}]})!.notes;
    expect(notes.length, 2, 'a note without a day is dropped');
    expect(notes[1].kind, 'Observation', 'an unknown kind reads as an observation');
  });

  test('the calendar marks each arm dosed by a dot, and Start for an arm that starts late', async () => {
    const g = (group: string, schedule: string) => ({group, compound: '', dose: null, unit: '', route: '', schedule});
    const owner = Scope.root('calendar test');
    const grid = (groups: ReturnType<typeof g>[]) => owner.build(() => weeksGrid({meta: {...emptyMeta(), groups},
      startWeekday: 0, last: 10, day: (_d, study) => ({marks: study}), onMenu: () => {}, onNote: () => {}}));
    const cellOf = (el: HTMLElement, day: number) => [...el.querySelectorAll('.pmx-calendar-day')]
      .find((c) => c.firstElementChild?.textContent?.startsWith(`Day ${day}`)) as HTMLElement;
    const dots = (c: HTMLElement) => c.querySelectorAll('.pmx-calendar-dot').length;
    const has = (c: HTMLElement, t: string) => [...c.querySelectorAll('span')].some((s) => s.textContent === t);
    const both = grid([g('Control', ''), g('A', 'QD'), g('B', 'Q3D')]);
    expect(dots(cellOf(both, 1)), 1, 'QD alone on day 1');
    expect(dots(cellOf(both, 3)), 2, 'QD and Q3D on day 3');
    expect(has(cellOf(both, 3), 'Dose'), false, 'no Dose word');
    expect(has(cellOf(both, 3), 'Start'), false, 'all start on day 0: no Start');
    expect(cellOf(both, 3).getAttribute('aria-label')!.includes('Dose: A, B'), true, 'the text names the arms');
    const late = grid([g('A', 'QD'), g('B', 'QD from day 7')]);
    expect(has(cellOf(late, 7), 'Start'), true, 'B starts on day 7');
    expect(has(cellOf(late, 8), 'Start'), false);
    expect(dots(cellOf(late, 6)), 1);
    expect(dots(cellOf(late, 7)), 2);
    owner.dispose();
  });

  test('a day clicked: its context panel offers what its right click does (notes, doses, its date)', async () => {
    const group = {group: 'A', compound: '', dose: null, unit: '', route: '', schedule: 'QD'};
    let meta = {...emptyMeta(), groups: [group], start: '2026-03-02',
      notes: [{day: 3, kind: 'Observation', text: 'ulcer', animal: '', who: 'me'}]} as StudyMeta;
    const days = dayActions(() => meta, (m) => meta = m, () => [], {last: () => 10, logGiven: true});
    const before = grok.shell.o;
    days.select(3);
    await delay(50);
    try {
      // (a day is its own object: one element would not replace another as the current object)
      const panel = (grok.shell.o as PmxDay).details();
      const text = panel.textContent ?? '';
      for (const s of ['Add a note…', 'Remove note: Observation: ulcer', 'Skip a dose…', 'Log doses given…',
        'Change the date of day 3…', 'Every other day of the study follows from it', 'Observation: ulcer'])
        expect(text.includes(s), true, s);
      // what happens on the day first; copying the date in the right-click menu only
      expect(panel.querySelector('.pmx-object-title')?.textContent?.startsWith('Day 3'), true);
      expect(text.includes('Dose') && text.includes('Notes'), true, text);
      expect(text.includes('Copy the date'), false);
      // its links run the same actions: the note removed
      [...panel.querySelectorAll<HTMLElement>('a, .ui-link')]
        .find((a) => a.textContent?.startsWith('Remove note'))!.click();
      expect(meta.notes.length, 0);
    } finally {
      grok.shell.o = before;
    }
  });

  test('one reader of schedules: rates, counts, first days, on and off', async () => {
    const days = (s: string, last: number) => dosedDays(readSchedule(s), last).join(',');
    expect(days('QD', 6), '0,1,2,3,4,5,6', 'bare QD: every day to the last');
    expect(readSchedule('BID')!.perDay, 2);
    expect(days('QOD', 6), '0,2,4,6');
    expect(days('Q2D', 6), '0,2,4,6');
    expect(days('Q3D', 9), '0,3,6,9');
    expect(days('Q4D', 8), '0,4,8');
    expect(days('QW', 14), '0,7,14');
    expect(days('BIW', 13), '0,3,7,10');
    expect(days('Q2W', 28), '0,14,28');
    expect(days('QD×5', 13), '0,1,2,3,4');
    expect(days('qdx3', 13), '0,1,2', 'a count with no space');
    expect(days('QD from day 3', 6), '3,4,5,6');
    expect(days('Q3D d3–d16', 30), '3,6,9,12,15');
    expect(days('days 0–4', 13), '0,1,2,3,4', 'days alone: daily');
    expect(days('5 on / 2 off', 13), '0,1,2,3,4,7,8,9,10,11');
    expect(days('5 on/2 off ×2', 30), '0,1,2,3,4,7,8,9,10,11', 'two cycles');
    expect(readSchedule('x3'), null, 'a count alone has no rate');
    expect(readSchedule('osmotic pump')!.kind, 'continuous');
    expect(readSchedule('', 'in drinking water')!.kind, 'continuous', 'from the route');
    expect(readSchedule('Q1 w/ PD-1'), null, 'the unit right after the number');
    // never continuous from a group's name: a diet in the name is not a dose
    const g0 = (group: string, schedule = '', route = '') => ({group, compound: '', dose: null, unit: '', route,
      schedule});
    expect([...dosingDays([g0('High-fat diet + Cpd-A 10 mg/kg Q3D')], 6).keys()].join(','), '0,3,6');
    expect(dosingDays([g0('Diet control'), g0('Vehicle (chow)'), g0('Cpd-A in drinking water')], 6).size, 0);
    expect(dosingDays([g0('Cpd-A', '', 'drinking water')], 6).size, 7, 'continuous when the route says so');
    // a parenthesis that is the schedule's days stays with the schedule
    const par = parseArmName('Cpd-A 10 mg/kg QD (d0-d4)')[0];
    expect(par.key, 'cpda');
    expect(dosedDays(readSchedule(par.schedule), 13).join(','), '0,1,2,3,4');
    // the arm's name reads the same schedules
    const a = parseArmName('Cpd-A 30 mg/kg QD from day 3 p.o.')[0];
    expect([a.drug, a.schedule, a.route].join('|'), 'Cpd-A|QD from day 3|p.o.');
    expect(parseArmName('Cpd-B 10 mg/kg i.p. Q3D×4, d3–d16')[0].schedule, 'Q3D×4, d3–d16');
    expect(parseArmName('Cpd-C 25 mg/kg p.o. 5 on / 2 off')[0].drug, 'Cpd-C');
    // a count after the rate stays with the schedule, never in the drug
    for (const n of ['Cpd-A 30 mg/kg p.o. QD×14', 'Cpd-A 30 mg/kg p.o. QW×3', 'Cpd-A 30 mg/kg p.o. BID x 14'])
      expect(parseArmName(n)[0].key, 'cpda', n);
    // the exposure's doses a day and the calendar's days agree (two weeks: QD 14, BID 28, BIW 4)
    for (const [sc, n] of [['QD', 14], ['BID', 28], ['BIW', 4]] as [string, number][]) {
      const s = readSchedule(sc)!;
      expect(dosedDays(s, 13).length * s.perDay, n, sc);
      expect(Math.abs(dosesPerDay(sc)! * 14 - n) < 1e-9, true, `${sc}: doses a day × 14 days`);
    }
    // the calendar reads them too
    const g = (group: string, schedule: string) => ({group, compound: '', dose: null, unit: '', route: '', schedule});
    expect(dosingDays([g('A', 'QD')], 3).size, 4, 'a bare QD marks every day');
  });

  test('no dosing at weekends: skipped, moved, kept', async () => {
    const monday = {start: '2026-11-09', weekendsOff: true};
    const days = (s: string, last: number | null, rules: {start: string | null, weekendsOff: boolean}) =>
      dosedDays(readSchedule(s), last, rules);
    expect(days('QD', 13, monday).length, 10, 'QD over two weeks: 10 doses');
    expect(days('QD', 13, {start: null, weekendsOff: true}).length, 14, 'no date of day 0: no effect');
    expect(days('QD', 13, {...monday, weekendsOff: false}).length, 14);
    expect(days('Q3D×4', null, monday).join(','), '0,3,7,9', 'Sunday to Monday, the count kept');
    expect(days('QW×3', null, {start: '2026-11-14', weekendsOff: true}).join(','), '2,6,13',
      'day 0 a Saturday: not before the first day, so Monday; later Saturdays to Friday');
    // Friday's dose taken: Saturday's goes the other way (never two doses of one arm on one day)
    expect(days('4 times a week ×4', null, {start: '2026-11-13', weekendsOff: true}).join(','), '0,3,4,5');
    expect(days('osmotic pump', 13, monday).length, 14, 'continuous exposure untouched');
    expect(days('5 on/2 off ×2', null, {start: '2026-11-11', weekendsOff: true}).length, 10,
      'on and off from a Wednesday: its weekdays');
    expect(dosesPerDay('QD', monday), 5 / 7);
    expect(dosesPerDay('Q3D', monday), 1 / 3, 'a dose moved is still a dose');
    expect(dosesPerDay('QD', {start: null, weekendsOff: true}), 1);
    // ×n counts doses: they go on past the weekend, whatever day 0's weekday
    const friday = {start: '2026-11-13', weekendsOff: true};
    expect(days('QD×3', null, friday).join(','), '0,3,4', 'Friday, Monday, Tuesday');
    expect(dosesPerDay('QD×3', friday), dosesPerDay('QD×3', monday), 'the same daily rate from any weekday');
    // a moved dose stays inside the study and the schedule's days, as the rate counts it
    expect(days('Q3D', 6, monday).join(','), '0,3,4', 'Sunday on the last day: Monday is past it, so Friday');
    expect(days('Q3D d0-d6', null, monday).join(','), '0,3,4', 'and past the schedule\'s last day too');
    expect(days('Q3D d3-d15', null, monday).join(','), '3,7,9,11,15', 'Sunday 6 → Monday 7, Saturday 12 → Friday 11');
    expect(days('QD', 13, {start: '2026-02-31', weekendsOff: true}).length, 14, 'a date that does not exist: no date');
    expect(parseMeta({start: '2026-02-31'})!.start, null);
    expect(parseMeta({weekendsOff: true})!.weekendsOff, true);
    expect(emptyMeta().weekendsOff, false);
  });

  test('dose exceptions: skipped and added, one shape for every reader', async () => {
    const g = (group: string, schedule: string) => ({group, compound: '', dose: null, unit: '', route: '', schedule});
    let m: StudyMeta = {...emptyMeta(), groups: [g('A', 'QD×5'), g('B', 'Q3D')]};
    m = withException(m, 2, 'A', 'skip', 'unwell');
    m = withException(m, 4, 'B', 'add');
    const days = dosingDays(m.groups, 6, m);
    expect((days.get(2) ?? []).includes('A'), false, 'a skipped dose is not dosed');
    expect(days.get(4)!.join(','), 'A,B', 'an added dose is');
    expect(m.groups[0].schedule, 'QD×5', 'the schedule as written never changes');
    expect(excepted([0, 1, 2], 'A', {start: null, exceptions: [{day: 9, group: 'A', change: 'add', note: '',
      who: '', when: ''}]}, 6).join(','), '0,1,2', 'an added day past the last is not dosed');
    // a removed arm leaves no doses behind (designer N1): its details and exceptions go with it
    const removed = onlyGroups(m, ['A']);
    expect([...dosingDays(removed.groups, 6, removed).values()].flat().includes('B'), false, 'B removed: not dosed');
    expect(removed.exceptions.some((e) => e.group === 'B'), false);
    expect(onlyGroups(m, ['A', 'B']) === m, true, 'nothing removed: the same details');
    const unread = withException({...emptyMeta(), groups: [g('C', 'whenever')]}, 1, 'C', 'add');
    expect(dosingDays(unread.groups, 6, unread).size, 0, 'no schedule read: an added dose does not dose it');
    expect(armDoses('Cpd-A 20 mg/kg QD×5', withException(emptyMeta(), 9, 'Cpd-A 20 mg/kg QD×5', 'skip'), 10, 10,
      0.025)!.text.includes('skipped'), false, 'a skip on a day not dosed changes nothing and is not counted');
    expect(withException(m, 2, 'A', null).exceptions.length, 1, 'undone: the other stays');
    expect(withException(m, 2, 'A', 'add').exceptions.filter((e) => e.day === 2).length, 1, 'one a day and group');
    // a skipped dose lowers the total, nothing is redistributed
    expect(armDoses('Cpd-A 20 mg/kg QD×5', {...emptyMeta(), exceptions: [{day: 1, group: 'Cpd-A 20 mg/kg QD×5',
      change: 'skip', note: '', who: '', when: ''}]}, 10, 10, 0.025)!.text.startsWith('4 doses · 80 mg/kg'), true);
    expect(doseForTotal('Cpd-A 20 mg/kg QD×5', emptyMeta(), 10, 100), 20, 'the dose that gives 100 mg/kg in 5');
    expect(doseForTotal('Vehicle', emptyMeta(), 10, 100), null);
    // doses given, as Monitor logs them
    const logged = withGiven(withGiven(m, 0, ['A', 'B']), 0, ['A']);
    expect(logged.given.map((x) => x.group).join(','), 'A', 'the day\'s log is what was ticked');
    const arm = 'Cpd-A 20 mg/kg QD×5';
    expect(armDoses(arm, withGiven({...emptyMeta()}, 0, [arm]), 10, 10, 0.025, 2)!.text
      .endsWith('1 of 3 dosing days so far logged as given'), true);
    // read back from a file, and renamed with their group
    const back = parseMeta(JSON.parse(JSON.stringify(logged)))!;
    expect(back.exceptions.length, 2);
    expect(back.given.length, 1);
    expect(parseMeta({exceptions: [{day: 1.5, group: 'A', change: 'skip'}, {day: 1, group: 'A', change: 'x'}]})!
      .exceptions.length, 0, 'not a day or not a change: dropped');
    const renamed = renamedDetails(m, new Map([['A', 'A2']]));
    expect(renamed.exceptions.find((e) => e.day === 2)!.group, 'A2');
  });

  test('what an arm adds up to: doses, dose in all, compound', async () => {
    const m = emptyMeta();
    expect(armDoses('Cpd-A 20 mg/kg QD×14', m, 30, 10, 0.025)!.text,
      '14 doses · 280 mg/kg in total · about 84 mg of compound for the group', '70 mg and 20% for losses');
    expect(armDoses('Cpd-A 20 mg/kg QD', {...m, start: '2026-11-09', weekendsOff: true}, 13, 10, null)!.text,
      '10 doses · 200 mg/kg in total · about 60 mg of compound for the group, assuming 25 g mice');
    expect(armDoses('Cpd-A 20 umol/kg QD×14', m, 30, 10, 0.025)!.text, '14 doses', 'not in mg/kg: the count alone');
    expect(armDoses('Cpd-A QD×14', m, 30, 10, 0.025)!.text, '14 doses', 'no dose: the count alone');
    expect(armDoses('Cpd-A 20 mg/kg QD×14', m, 30, 10, 0.025, 5)!.text,
      '6 of 14 doses so far · 120 of 280 mg/kg · about 84 mg of compound for the group');
    expect(armDoses('Vehicle', m, 30, 10, 0.025), null);
    expect(armDoses('Cpd-A 20 mg/kg QD×14', {...m, start: '2026-11-09', weekendsOff: true}, 30, 10, 0.025)!.text
      .startsWith('14 doses'), true, 'a count of doses stays at weekends off');
    // the animals still in the study and their last weights: gone before the group's last day (dead or removed)
    // and excluded rows left out; grams by the column's unit, else by a mean of 10 or more
    const row = (animal: string, day: number, weight: number | null, excluded = false) =>
      ({group: 'A', animal, day, weight, excluded});
    const rows = [row('M1', 0, 20), row('M1', 7, 22), row('M2', 0, 21), row('M3', 0, 24), row('M3', 7, 30, true)];
    const a = groupWeights(rows, null).get('A')!;
    expect(a.animals, 1);
    expect(Math.abs(a.kg! - 0.022) < 1e-12, true, String(a.kg));
    expect(groupWeights(rows.map((r) => ({...r, weight: r.weight && r.weight / 1000})), null).get('A')!.kg! < 0.03,
      true, 'under 10: kilograms');
    expect(Math.abs(groupWeights(rows, 'kg').get('A')!.kg! - 22) < 1e-9, true, 'the column says kg');
  });

  test('a saved study keeps an empty day empty, not day 0', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M1']),
      DG.Column.fromStrings('group', ['A', 'A']), DG.Column.fromList('double', 'day', [3, null]),
      DG.Column.fromList('double', 'volume', [100, null])]);
    const back = tableOf(parseStudy(JSON.stringify(studyOf(df, null))));
    expect(back.get('day', 0), 3);
    expect(back.col('day')!.isNone(1), true, 'the empty day came back as a value');
    expect(back.col('volume')!.isNone(1), true, 'the empty volume came back as a value');
  });
  test('a study opened from a table is listed by its groups, unless they are blinded', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M2', 'M3']),
      DG.Column.fromStrings('group', ['Vehicle', 'Cpd-A 10 mg/kg', 'Vehicle']), DG.Column.fromList('double', 'day',
        [0, 0, 0]), DG.Column.fromList('double', 'volume', [100, 90, 95])]);
    const entry = () => entryOf(JSON.stringify(studyOf(df, null)), 1)!;
    expect([entry().compounds, entry().treatment].join('|'), 'Vehicle, Cpd-A 10 mg/kg|Cpd-A 10 mg/kg',
      'its own groups');
    setMetaOn(df, {...emptyMeta(), blinded: true});
    expect([entry().compounds, entry().treatment].join('|'), '|', 'blinded: no group named');
    addLog(df, REVEALED);
    expect(entry().treatment, 'Cpd-A 10 mg/kg', 'revealed: named again');
    setMetaOn(df, {...emptyMeta(), groups: [{group: 'Cpd-A 10 mg/kg', compound: 'Cpd-A', dose: 10, unit: 'mg/kg',
      route: '', schedule: ''}]});
    expect(entry().compounds, 'Cpd-A', 'the details first, the groups only without them');
  });
  test('an arm\'s name read into its details, shown, kept and listed', async () => {
    // the line under the name: one drug, a combination, a control
    expect(parseLine('Cpd-A 30 mg/kg QD PO'), 'Cpd-A \u00b7 30 mg/kg \u00b7 QD \u00b7 p.o.');
    expect(parseLine('Cpd-A 30 mg/kg QD PO + anti-PD-1 10 mg/kg i.p. Q3D'),
      'Cpd-A \u00b7 30 mg/kg \u00b7 QD \u00b7 p.o. + anti-PD-1 \u00b7 10 mg/kg \u00b7 Q3D \u00b7 i.p.');
    expect(parseLine('Vehicle'), 'control');
    // the details' empty cells from the name, a typed one kept; a combination one row per drug
    const name = 'Cpd-A 30 mg/kg QD PO';
    const typed = {...emptyMeta(), groups: [{group: name, compound: '', dose: null, unit: '', route: 'i.v.',
      schedule: ''}]};
    const [g] = withGroups(typed, [name]).groups;
    expect(JSON.stringify([g.compound, g.dose, g.unit, g.route, g.schedule]), '["Cpd-A",30,"mg/kg","i.v.","QD"]');
    expect(parseLine(name, withGroups(typed, [name])), 'Cpd-A \u00b7 30 mg/kg \u00b7 QD \u00b7 i.v. (route from the ' +
      'group details)');
    const combo = 'Cpd-A 30 mg/kg + anti-PD-1 10 mg/kg';
    expect(withGroups(emptyMeta(), [combo]).groups.map((g) => g.compound).join(), 'Cpd-A,anti-PD-1');
    expect(parseLine(combo, withGroups(emptyMeta(), [combo])).includes(' + '), true);
    // a renamed group: the details typed under its old name stay in the meta; a typed unit is kept
    const renamed = withGroups(typed, ['Cpd-A 50 mg/kg QD PO']).groups;
    expect(renamed.find((r) => r.group === name)?.route, 'i.v.');
    expect(fromName({group: 'Cpd-A 30 mg/m2', compound: '', dose: null, unit: 'mg/kg', route: '', schedule: ''})
      .unit, 'mg/kg');
    // saved and reopened: the details travel with the study, and are read first
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M2']),
      DG.Column.fromStrings('group', ['Vehicle', 'Low']), DG.Column.fromList('double', 'day', [0, 0]),
      DG.Column.fromList('double', 'volume', [100, 90])]);
    setMetaOn(df, {...emptyMeta(), groups: [{group: 'Low', compound: 'Cpd-A', dose: 10, unit: 'mg/kg', route: 'p.o.',
      schedule: 'QD'}]});
    const back = metaOf(tableOf(parseStudy(JSON.stringify(studyOf(df, null)))));
    expect(parseLine('Low', back), 'Cpd-A · 10 mg/kg · QD · p.o. (compound, dose, route, ' +
      'schedule from the group details)');
    // an older study, no details: read from the name as before
    expect(JSON.stringify(armParts(name, null)), JSON.stringify(parseArmName(name)));
    expect(JSON.stringify(armParts(name, emptyMeta())), JSON.stringify(parseArmName(name)));
    // the Overview: the drugs from the details, else the names; found by them
    const entry = entryOf(JSON.stringify(studyOf(df, null)), 1)!;
    expect(entry.drugs, 'Cpd-A');
    expect(matches(entry, 'cpd-a'), true);
    const old = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1']),
      DG.Column.fromStrings('group', ['Cpd-B 5 mg/kg']), DG.Column.fromList('double', 'day', [0]),
      DG.Column.fromList('double', 'volume', [100])]);
    expect(entryOf(JSON.stringify(studyOf(old, null)), 1)!.drugs, 'Cpd-B');
  });
  test('a study keeps its tumor model without a plan; with one, the plan\u2019s model is read', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1']), DG.Column.fromStrings('group', ['A']),
      DG.Column.fromList('double', 'day', [0]), DG.Column.fromList('double', 'volume', [100])]);
    // a file saved before the field (none): read as before, no model of its own
    const old = studyOf(df, null);
    expect('line' in old, false);
    expect([lineOf(tableOf(parseStudy(JSON.stringify(old)))), entryOf(JSON.stringify(old), 1)!.line].join('|'),
      '|', 'an older file: no model');
    // saved now, without a plan: the model chosen comes back on the table and in the Overview
    df.setTag(LINE_TAG, 'B16F10');
    const text = JSON.stringify(studyOf(df, null));
    expect(lineOf(tableOf(parseStudy(text))), 'B16F10');
    expect(entryOf(text, 1)!.line, 'B16F10');
    // a field a version does not know is kept, not refused (as an older PMX Studio reads this one); a model that is
    // not a name is dropped
    expect((parseStudy(JSON.stringify({...JSON.parse(text), future: 1})) as any).future, 1);
    expect(parseStudy(JSON.stringify({...JSON.parse(text), line: 7})).line === undefined, true);
    // with a plan of another model, the plan's is read: it is what the study was sized and forecast for
    const study = await exampleStudy();
    const planLine = parsePlan(study.getTag(PLAN_TAG)!).predicted.line;
    study.setTag(LINE_TAG, planLine === 'B16F10' ? 'CT26' : 'B16F10');
    const both = JSON.stringify(studyOf(study, null));
    expect([lineOf(tableOf(parseStudy(both))), entryOf(both, 1)!.line].join('|'), `${planLine}|${planLine}`);
  });
  test('a study saved from Analyze is listed as finished, and stays so when saved again from its table', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1']), DG.Column.fromStrings('group', ['A']),
      DG.Column.fromList('double', 'day', [0]), DG.Column.fromList('double', 'volume', [100])]);
    expect(entryOf(JSON.stringify(studyOf(df, null)), 1)!.stage, 2, 'followed');
    df.setTag(FINISHED_TAG, 'true');
    const text = JSON.stringify(studyOf(df, null));
    expect(entryOf(text, 1)!.stage, 3, 'finished');
    expect(entryOf(JSON.stringify(studyOf(tableOf(parseStudy(text)), null)), 1)!.stage, 3, 'reopened and saved');
  });
  test('a study saved from Analyze keeps its method, endpoint and evaluation day (QA 0.1.3 r4: Holm came back as ' +
    'Bonferroni)', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1']), DG.Column.fromStrings('group', ['A']),
      DG.Column.fromList('double', 'day', [0]), DG.Column.fromList('double', 'volume', [100])]);
    expect(studyOf(df, null).analysis === undefined, true, 'none set: none in the file');
    setAnalysisOn(df, {correction: 'holm', endpoint: 1200, evaluationDay: 14});
    const text = JSON.stringify(studyOf(df, null));
    expect(JSON.stringify(JSON.parse(text).analysis), '{"correction":"holm","endpoint":1200,"evaluationDay":14}');
    const back = tableOf(parseStudy(text));
    expect(JSON.stringify(analysisFieldsOf(back)), '{"correction":"holm","endpoint":1200,"evaluationDay":14}',
      'reopened: the form gets them back');
    // a file edited by hand: only what may be there
    const edited = JSON.stringify({...JSON.parse(text), analysis: {correction: 'magic', endpoint: -1,
      evaluationDay: 7}});
    expect(JSON.stringify(analysisFieldsOf(tableOf(parseStudy(edited)))), '{"evaluationDay":7}');
  });
  test('a followed study goes by its plan\'s name; a plan with no name leaves the table\'s', async () => {
    const study = await exampleStudy();
    const plan = JSON.parse(study.getTag(PLAN_TAG)!);
    study.setTag(PLAN_TAG, JSON.stringify({...plan, name: 'MC38 full study'}));
    expect(studyOf(study, null).name, 'MC38 full study');
    study.setTag(PLAN_TAG, JSON.stringify({...plan, name: ''}));
    expect(studyOf(study, null).name, study.name);
  });

  /** A studies folder in memory: name -> the ID saved there (null: a file with none). */
  const folder = (files: {[name: string]: string | null}) => ({
    owner: async (name: string) => name in files ? files[name] : undefined,
    save: async (names: string[], id: string) => {
      const name = await freeName(names, id, async (n) => n in files ? files[n] : undefined);
      files[name] = id;
      return name;
    },
  });

  test('a study name taken by another study is never replaced: the table\'s name, then numbered', async () => {
    const {owner} = folder({'P': 'A', 'T': 'A', 'T (2)': 'C', 'X': null});
    expect(await freeName(['P', 'T'], 'B', owner), 'T (3)', 'plan and table names taken: numbered past the copies');
    expect(await freeName(['P', 'U'], 'B', owner), 'U', 'the plan\'s taken: the table\'s');
    expect(await freeName(['X'], 'B', owner), 'X', 'a file with no ID is no one\'s: replaced');
    expect(await freeName(['P'], 'A', owner), 'P', 'its own file: replaced');
    expect(await freeName(['', ' '], 'B', owner), 'Study', 'no name at all');
  });

  test('a copy past a deleted numbered file is found: the study is not saved twice', async () => {
    const {owner} = folder({'T': 'A', 'T (3)': 'B'});
    expect(await freeName(['T'], 'B', owner), 'T (3)', 'its own file, past the gap');
    expect(await freeName(['T'], 'C', owner), 'T (2)', 'a new study fills the gap');
  });

  test('two studies followed from one plan are saved apart, each again in its own file', async () => {
    const f = folder({});
    expect(await f.save(['MC38 plan', 'mc38'], 'A'), 'MC38 plan');
    expect(await f.save(['MC38 plan', 'mc38'], 'B'), 'mc38');
    expect(await f.save(['MC38 plan', 'mc38'], 'C'), 'mc38 (2)');
    expect([await f.save(['MC38 plan', 'mc38'], 'B'), await f.save(['MC38 plan', 'mc38'], 'C')].join('|'),
      'mc38|mc38 (2)', 'saved again: the same file, not a new copy');
  });

  test('a study saved under its table\'s name opens as before and is saved again in that file', async () => {
    const study = await exampleStudy();
    const plan = JSON.parse(study.getTag(PLAN_TAG)!);
    study.setTag(PLAN_TAG, JSON.stringify({...plan, name: 'MC38 full study'}));
    setMetaOn(study, {...emptyMeta(), id: 'PMX-2026-001'});
    // as an earlier PMX Studio saved it: by its table's name
    const old = parseStudy(JSON.stringify({...studyOf(study, null), name: 'mc38-day6'}));
    const back = tableOf(old);
    expect(back.name, 'mc38-day6');
    expect(back.rowCount, study.rowCount);
    const {owner} = folder({'mc38-day6': 'PMX-2026-001'});
    expect(await freeName([studyOf(back, null).name, back.name], 'PMX-2026-001', owner), 'mc38-day6',
      'not a second file under the plan\'s name');
  });

  test('a blinded study names no arm or compound outside itself until revealed', async () => {
    const study = await exampleStudy();
    const plan = parsePlan(study.getTag(PLAN_TAG)!);
    setMetaOn(study, {...emptyMeta(), id: 'PMX-2026-099', blinded: true, groups: plan.predicted.arms.map((a) =>
      ({group: a.name, compound: a.name.split(' ')[0], dose: null, unit: '', route: '', schedule: ''}))});
    // the arms and their compounds (the study's own name is the user's, and stays)
    const names = [...plan.predicted.arms.map((a) => a.name), ...plan.predicted.arms.map((a) => a.name.split(' ')[0])]
      .filter((n) => n);
    const said = () => {
      const text = JSON.stringify(studyOf(study, null));
      const e = entryOf(text, 1)!;
      // the Overview's row, its search, its index row and the file's preview
      return {e, text: [e.compounds, e.treatment, e.line, ...e.tags, ...studySummary(parseStudy(text))].join(' | ')};
    };
    const hidden = said();
    const leaks = names.filter((n) => hidden.text.includes(n) || matches(hidden.e, n));
    expect(leaks.length, 0, `blinded, named: ${leaks.join(', ')} in ${hidden.text}`);
    expect(matches(hidden.e, 'PMX-2026-099'), true, 'found by its ID');
    addLog(study, REVEALED);
    const compound = names[names.length - 1];
    expect(said().e.compounds.includes(compound) && matches(said().e, compound), true, 'revealed: named again');
    // a file whose blinding is in its plan's details only (no details of its own): blinded all the same
    const x = JSON.parse(JSON.stringify(studyOf(study, null)));
    delete x.meta;
    x.log = [];
    x.plan.study = {...emptyMeta(), blinded: true};
    const planOnly = entryOf(JSON.stringify(x), 1)!;
    const shownThere = [planOnly.compounds, planOnly.treatment, ...studySummary(parseStudy(JSON.stringify(x)))]
      .join(' | ');
    expect(names.filter((n) => shownThere.includes(n) || matches(planOnly, n)).join(', '), '',
      `blinded by its plan, named in ${shownThere}`);
  });
  test('a saved study keeps its marks and its details', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['M1', 'M2']),
      DG.Column.fromStrings('group', ['A', 'B']), DG.Column.fromList('double', 'day', [0, 0]),
      DG.Column.fromList('double', 'volume', [100, 90])]);
    mark(df, [1], 'excluded', 'misread', 'M2 on day 0');
    const meta = {...emptyMeta(), id: 'PMX-2026-007', start: '2026-11-02', tags: ['IO'],
      groups: [{group: 'B', compound: 'anti-PD-1', dose: 10, unit: 'mg/kg', route: 'i.p.', schedule: 'Q3D×4'}]};
    setMetaOn(df, meta);
    const back = tableOf(parseStudy(JSON.stringify(studyOf(df, null))));
    expect(back.col(EXCLUDED)!.toList().join(','), 'false,true');
    expect(back.get(MARK_REASON, 1), 'misread');
    expect(withoutExcluded({study: back}).study.rowCount, 1);
    expect(JSON.stringify(metaOf(back)), JSON.stringify(meta));
    // the library: the study and a plan with the same ID are one entry, at the study's stage; found by a compound
    const study = entryOf(JSON.stringify(studyOf(df, null)), 2)!;
    const plan = {...study, kind: 'plan' as const, stage: 1 as const, updated: 3};
    const listed = libraryOf([plan, study, null]);
    expect(listed.length, 1);
    expect(listed[0].kind, 'study');
    expect(matches(listed[0], 'pd-1 pmx-2026'), true);
    expect(matches(listed[0], 'cisplatin'), false);
    expect(study.treatment, 'anti-PD-1 10 mg/kg i.p. Q3D×4', 'each compound with its dose, route and schedule');
    // a combination's drugs joined by " + ", the groups by "; " (a " + " between groups reads as a combination)
    const arm = (group: string, compound: string) => ({group, compound, dose: null, unit: '', route: '',
      schedule: ''});
    expect(treatmentOf({...emptyMeta(), groups: [arm('A', 'Cpd-A'), arm('AB', 'Cpd-A'), arm('AB', 'anti-PD-1'),
      arm('B', 'anti-PD-1')]}), 'Cpd-A; Cpd-A + anti-PD-1; anti-PD-1');
    // rows without a group's name are groups of their own, not one combination (QA 0.1.3 #8)
    expect(treatmentOf({...emptyMeta(), groups: [arm('', 'Cpd-A'), arm('', 'Cpd-B')]}), 'Cpd-A; Cpd-B');
    expect(compoundsOf({...emptyMeta(), groups: [arm('A', 'Trametinib'), arm('B', 'trametinib')]}), 'Trametinib');
    expect([study.start, study.day, study.end].join(), '2026-11-02,0,', 'day 0, the last day measured, no plan');
    expect(parseMeta({id: 3, start: '2 Nov', groups: [{group: 'A', dose: 'ten'}]})!.groups[0].dose, null,
      'a damaged field is dropped, not trusted');
    // blinded entry: the groups hidden until revealed; each mark says whether its author knew them
    expect(groupsShown(df), true, 'not blinded: the groups are shown');
    setMetaOn(df, {...meta, blinded: true});
    expect(groupsShown(df), false);
    mark(df, [0], 'flagged', 'odd shape', 'M1 on day 0');
    addLog(df, REVEALED);
    expect(groupsShown(df), true, 'revealed, for good');
    mark(df, [0], null, '', 'M1 on day 0');
    const said = logOf(df).map((e) => e.what);
    expect(said.some((w) => w.startsWith('Flagged M1') && w.endsWith('(blinded)')), true);
    expect(said[said.length - 1].endsWith('(groups revealed)'), true);
    expect(metaOf(tableOf(parseStudy(JSON.stringify(studyOf(df, null)))))!.blinded, true);
  });

  test('an animal picked in one chart is selected in the study, and in every table of animals', async () => {
    const s = study();
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', s.animal as string[]),
      DG.Column.fromStrings('group', s.group as string[]), DG.Column.fromList('double', 'day', s.day),
      DG.Column.fromList('double', 'volume', s.volume)]);
    const sel = new AnimalSelection();
    sel.attach({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null});
    const settled = () => new Promise((r) => setTimeout(r, 50));
    const animals = () => [...sel.animals.peek()].sort().join(',');
    sel.pick('M2');
    await settled();
    expect(df.selection.trueCount, 3, 'its three rows');
    expect(animals(), 'M2');
    sel.pick('M4', new MouseEvent('click', {ctrlKey: true}));
    await settled();
    expect(animals(), 'M2,M4', 'Ctrl adds');
    // a table of animals follows, and its own selection reaches the study
    const each = DG.DataFrame.fromColumns([DG.Column.fromStrings('Animal', ['M1', 'M2', 'M3', 'M4'])]);
    const stop = mirrorSelection(each, 'Animal', sel);
    expect(each.selection.trueCount, 2);
    each.selection.init((i) => i === 0);
    await settled();
    expect(animals(), 'M1', 'a selection on the derived table');
    stop();
    sel.pick('M1');
    await settled();
    expect(df.selection.trueCount, 0, 'picked again alone: none');
    sel.dispose();
  });

  test('"values to check" selects their rows in the table and makes the first current', async () => {
    const s = study();
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', s.animal as string[]),
      DG.Column.fromStrings('group', s.group as string[]), DG.Column.fromList(DG.TYPE.FLOAT, 'day', s.day),
      DG.Column.fromList(DG.TYPE.FLOAT, 'volume', s.volume)]);
    const n = {study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null};
    const found = findingsOf(n);
    expect(found.size > 0, true, 'the example has a value to check');
    const line = studyDataLine(() => n);
    try {
      line.refresh();
      line.root.querySelector<HTMLElement>('.ui-link')!.click();
      expect(df.selection.trueCount, found.size);
      expect([...found.keys()].every((i) => df.selection.get(i)), true);
      expect(df.currentRowIdx, Math.min(...found.keys()));
    } finally {
      line.dispose();
      [...grok.shell.tableViews].find((v) => v.dataFrame.id === df.id)?.close();
    }
  });

  test('whole days show as whole numbers; a fractional one keeps its decimals', async () => {
    const day = DG.Column.fromList(DG.TYPE.FLOAT, 'day', [0, 3, 7]);
    dayFormat(day);
    expect(day.getTag(DG.TAGS.FORMAT), '0');
    day.set(1, 3.5);
    dayFormat(day);
    expect(!day.getTag(DG.TAGS.FORMAT), true, 'a fractional day: decimals back');
    const own = DG.Column.fromList(DG.TYPE.FLOAT, 'day', [0, 3]);
    own.meta.format = '0.0';
    dayFormat(own);
    expect(own.getTag(DG.TAGS.FORMAT), '0.0', 'a format of the user\'s own');
  });

  test('outlier rules out of range are refused with the allowed range, never saved', async () => {
    expect(rulesProblem(DEFAULT_RULES), null);
    const bad = {...DEFAULT_RULES, jump: {on: true, fold: 0.5}, group: {on: true, mads: 3.5, minAnimals: -5}};
    expect(rulesProblem(bad), 'Times larger or smaller: from 1.5 to 20; In groups of at least: a whole number from ' +
      '3 to 50.');
    expect(rulesProblem({...DEFAULT_RULES, weightDrop: {on: false, percent: NaN}}),
      'Loss since the last weighing (%): from 1 to 50.', 'an emptied threshold is not a default');
  });
});
