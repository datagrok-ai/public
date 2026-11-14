import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_RULES, findOutliers, StudyColumns} from '../outlier-rules';
import {AnimalSelection, animalMarks, appendPasted, caliperVolume, emptyStudy, exclusionsOf, mark, mirrorSelection,
  pastedNumber, withoutExcluded} from '../study-table';
import {EXCLUDED, FLAGGED, LINE_TAG, MARK_REASON, REVEALED, STUDY_TAG, WEIGHT_TAG, addLog, groupsShown, lineOf, logOf,
  freeName, parseStudy, roleColumns, studyOf, tableOf} from '../study-file';
import {PLAN_TAG, parsePlan} from '../plan-file';
import {exampleStudy} from './example-study';
import {filePreview, studySummary} from '../file-preview';
import {EXAMPLE, isBodyWeight} from '../ui';
import {_package} from '../package';
import {dateOf, dosingDays, emptyMeta, metaOf, parseMeta, setMetaOn, startFrom, weekdayOf} from '../study-meta';
import {entryOf, libraryOf, matches} from '../library';
import {spiderLines} from '../spider-chart';
import {dayToday, icsOf} from '../calendar';
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
    expect([...doses.keys()].sort((a, b) => a - b).join(','), '0,1,3,6,7,9,14');
    expect(doses.get(0)!.join(','), 'A,B,D', 'every group on day 0; BIW is not read');
    expect(doses.get(9)!.join(','), 'A');
    expect(dosingDays([g('A', 'Q7D×10')], 20).size, 3, 'nothing past the last day');
    expect(dateOf('2026-11-09', 22)!.toISOString().slice(0, 10), '2026-12-01');
    expect(dateOf(null, 3), null);
    expect(startFrom('2026-11-12', 3), '2026-11-09');
    expect(weekdayOf('2026-11-09'), 0, 'a Monday');
    expect(dayToday('2026-11-09', new Date(2026, 10, 21, 23, 30)), 12, 'today, by the local date');
    expect(dayToday('2026-11-09', new Date(2026, 10, 2)), -7, 'before day 0');
    expect(dayToday(null), null);
    const ics = icsOf('S', {...emptyMeta(), start: '2026-11-09', groups: [g('A', 'QD×2')]}, 3, [[0, 'Randomize']]);
    expect(ics.split('BEGIN:VEVENT').length - 1, 3, 'randomize and two doses, all-day');
    expect(ics.includes('DTSTART;VALUE=DATE:20261110'), true);
    const notes = parseMeta({notes: [{day: 3, kind: 'Incident', animal: 'M4', text: 'ulcerated', who: 'x', when: ''},
      {day: 'x', text: 'no day'}, {day: 5, kind: 'Party', text: 'odd kind'}]})!.notes;
    expect(notes.length, 2, 'a note without a day is dropped');
    expect(notes[1].kind, 'Observation', 'an unknown kind reads as an observation');
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
    expect([entry().compounds, entry().treatment].join('|'), 'Vehicle, Cpd-A 10 mg/kg|Vehicle, Cpd-A 10 mg/kg',
      'its own groups');
    setMetaOn(df, {...emptyMeta(), blinded: true});
    expect([entry().compounds, entry().treatment].join('|'), '|', 'blinded: no group named');
    addLog(df, REVEALED);
    expect(entry().treatment, 'Vehicle, Cpd-A 10 mg/kg', 'revealed: named again');
    setMetaOn(df, {...emptyMeta(), groups: [{group: 'Cpd-A 10 mg/kg', compound: 'Cpd-A', dose: 10, unit: 'mg/kg',
      route: '', schedule: ''}]});
    expect(entry().compounds, 'Cpd-A', 'the details first, the groups only without them');
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
});
