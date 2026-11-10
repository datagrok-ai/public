import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_RULES, findOutliers, StudyColumns} from '../outlier-rules';
import {AnimalSelection, animalMarks, appendPasted, caliperVolume, emptyStudy, exclusionsOf, mark, mirrorSelection,
  pastedNumber, withoutExcluded} from '../study-table';
import {EXCLUDED, FLAGGED, MARK_REASON, REVEALED, STUDY_TAG, addLog, groupsShown, logOf, parseStudy, studyOf, tableOf}
  from '../study-file';
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
