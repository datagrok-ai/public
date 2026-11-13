import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {FORMULAS, type SheetFrame, type TableRead, type VolumeFormula, cellsOfFrame, csvOf, dateOrder, deduceFormula,
  readStudy, roleOf, studySheets, studySummary,
  volumeSourceOf, waiting} from '../study-reader';
import {cellsOf} from '../wide-table';

const read = (csv: string, answers = {}) => readStudy(cellsOf(csv), answers);
const lines = (r: TableRead, blinded = false) => studySummary(r, blinded).join(' | ');
const near = (a: number | null, b: number) => a !== null && Math.abs(a - b) < 1e-6 * Math.max(1, Math.abs(b));

/** A long table with length and width (and height when `h`): two animals, two days. */
const caliper = (head: string, h = false, vol?: VolumeFormula) => {
  const rows: [string, string, number, number, number, number][] = [['m1', 'Vehicle', 0, 10, 8, 6],
    ['m1', 'Vehicle', 3, 12, 9, 7], ['m2', 'Cpd', 0, 11, 7, 5], ['m2', 'Cpd', 3, 11.5, 7.5, 6]];
  if (vol)
    rows.push(['m1', 'Vehicle', 1, 10.5, 8.5, 6.5], ['m1', 'Vehicle', 2, 11, 8.8, 6.8]);
  return head + '\n' + rows.map(([a, g, d, l, w, hh]) => [a, g, d, l, w, ...h ? [hh] : [],
    ...vol ? [Math.round(FORMULAS[vol].v(l, w, hh) * 100) / 100] : []].join(',')).join('\n') + '\n';
};

category('Study reader', () => {
  test('roles by name, closed lists, units from the name', async () => {
    expect(JSON.stringify(roleOf('Tumor Volume (mm3)')), JSON.stringify({role: 'volume', unit: 'mm3'}));
    expect(JSON.stringify(roleOf('Length_cm')), JSON.stringify({role: 'length', unit: 'cm'}));
    expect(JSON.stringify(roleOf('BW (kg)')), JSON.stringify({role: 'weight', unit: 'kg'}));
    expect(JSON.stringify(roleOf('Mouse ID')), JSON.stringify({role: 'animal', unit: null}));
    expect(roleOf('Tumor weight'), null);
    expect(roleOf('Body weight (% change)'), null);
    expect(roleOf('Volume (g)'), null);
  });

  test('a long table with length and width only: length × width² / 2, the others offered', async () => {
    const r = read(caliper('Animal,Group,Day,Length,Width'));
    expect(r.error, null);
    expect(waiting(r), false);
    expect(JSON.stringify(r.volume), JSON.stringify({formula: 'LW2/2', source: 'computed'}));
    expect(near(r.rows[0].volume, 10 * 8 * 8 / 2), true);
    expect(r.questions[0].options.map(([f]) => f).join(), 'LW2/2,pi/6*LW2');
    expect(lines(r), 'Read 2 groups: Vehicle 1 animal, Cpd 1 animal, days 0–3. | Volumes computed as length × ' +
      'width² / 2. | No unit in the column names: assumed mm. | Also kept: Length, Width.');
    const pi = read(caliper('Animal,Group,Day,Length,Width'), {formula: 'pi/6*LW2'});
    expect(near(pi.rows[0].volume, Math.PI / 6 * 640), true);
  });

  test('a long table with length, width and height: the formulas with a height offered too', async () => {
    const r = read(caliper('Animal,Group,Day,Length (mm),Width (mm),Height (mm)', true));
    expect(r.questions[0].options.map(([f]) => f).join(), 'LW2/2,pi/6*LW2,pi/6*LWH,LWH/2');
    const lwh = read(caliper('Animal,Group,Day,Length (mm),Width (mm),Height (mm)', true), {formula: 'LWH/2'});
    expect(near(lwh.rows[0].volume, 10 * 8 * 6 / 2), true);
    expect(lines(lwh).includes('assumed'), false);
  });

  test('length, width and the lab\'s volumes: its formula found, each of the four', async () => {
    for (const f of Object.keys(FORMULAS) as VolumeFormula[]) {
      const h = FORMULAS[f].height;
      const r = read(caliper(`Animal,Group,Day,Length,Width,${h ? 'Height,' : ''}Volume`, h, f));
      expect(JSON.stringify(r.volume), JSON.stringify({formula: f, source: 'deduced'}), f);
      expect(studySummary(r)[1], `Volumes match ${FORMULAS[f].label}.`, f);
    }
    expect(deduceFormula([{l: 10, w: 8, h: null, v: 999}]), null);
    // fewer than five rows: no formula is deduced, the volumes are as given
    expect(deduceFormula([1, 2, 3, 4].map((i) => ({l: 10 + i, w: 8, h: null, v: (10 + i) * 32}))), null);
    expect(JSON.stringify(volumeSourceOf({formula: 'LWH/2', source: 'deduced'})),
      JSON.stringify({formula: 'LWH/2', source: 'deduced'}));
    expect(volumeSourceOf({formula: 'L*W', source: 'computed'}) === undefined, true);
    expect(volumeSourceOf(undefined) === undefined, true);
  });

  test('volumes that follow no formula are used as given, and said so', async () => {
    const csv = caliper('Animal,Group,Day,Length,Width,Volume', false, 'LW2/2').replace(/,320\n/, ',400\n')
      .replace(/,486\n/, ',300\n').replace(/,171.5\n/, ',90\n');
    const r = read(csv);
    expect(JSON.stringify(r.volume), JSON.stringify({formula: null, source: 'given'}));
    expect(studySummary(r)[1], 'The volumes given do not follow a known formula from length and width: they are ' +
      'used as given.');
    expect(r.rows[0].volume, 400);
  });

  test('units from the names converted: cm to mm, cm³ and mL to mm³, kg to g', async () => {
    const cm = read(caliper('Animal,Group,Day,Length (cm),Width (cm)').replace(/,10,8\n/, ',1,0.8\n'));
    expect(near(cm.rows[0].volume, 320), true);
    expect(studySummary(cm).includes('Lengths converted from cm to mm.'), true, lines(cm));
    const ml = read('Animal,Group,Day,Volume (mL),BW (kg)\nm1,V,0,0.1,0.02\nm2,C,0,0.12,0.021\n');
    expect(ml.rows[0].volume, 100);
    expect(ml.rows[0].weight, 20);
    expect(lines(ml), 'Read 2 groups: V 1 animal, C 1 animal, days 0–0. | Volumes as given. | Volumes converted ' +
      'from mL to mm³. | Body weights converted from kg to g.');
  });

  test('no unit in the names: assumed and said; volumes that are no study\'s: the unit is asked', async () => {
    const r = read('Animal,Group,Day,Volume\nm1,V,0,100\nm2,C,0,120\n');
    expect(studySummary(r).includes('No unit in the column names: assumed mm³.'), true, lines(r));
    const small = 'Animal,Group,Day,Volume\nm1,V,0,0.1\nm2,C,0,0.12\n';
    const asked = read(small);
    expect(waiting(asked), true);
    expect(asked.rows.length, 0);
    expect(asked.questions[0].id, 'unit:volume');
    const answered = read(small, {'unit:volume': 'cm3'});
    expect(answered.rows[0].volume, 100);
    expect(lines(answered).includes('Volumes converted from cm³ to mm³.'), true, lines(answered));
    // lengths in cm with no unit: the volumes come out under 1 mm³, the lengths' unit is asked
    const cm = read('Animal,Group,Day,Length,Width\nm1,V,0,1,0.8\nm2,C,0,0.9,0.7\n');
    expect(cm.questions.map((q) => q.id).join(), 'formula,unit:length');
    expect(near(read('Animal,Group,Day,Length,Width\nm1,V,0,1,0.8\nm2,C,0,0.9,0.7\n', {'unit:length': 'cm'})
      .rows[0].volume, 320), true);
  });

  test('dates: ISO, day 0 the first measurement', async () => {
    const r = read('Animal,Group,Date,Volume (mm3)\nm1,V,2026-03-02,100\nm1,V,2026-03-05,150\nm2,C,2026-03-03,90\n');
    expect(r.rows.map((x) => x.day).join(), '0,3,1');
    expect(studySummary(r).includes('Day 0 is the first measurement, 2026-03-02.'), true, lines(r));
    // days and dates both: the days as given, the dates kept
    const both = read('Subject,Group,Date,Study Day,Volume (mm3)\nm1,V,2024-01-01,0,100\nm1,V,2024-01-04,3,150\n');
    expect(lines(both), 'Read 1 group: V 1 animal, days 0–3. | Volumes as given. | Also kept: Date.');
  });

  test('dates: dd/mm told by a day above 12, asked when no date tells', async () => {
    const dmy = read('Animal,Group,Date,Volume\nm1,V,28/02/2026,100\nm1,V,03/03/2026,150\n');
    expect(dmy.rows.map((x) => x.day).join(), '0,3');
    expect(dateOrder(['28/02/2026', '03/03/2026']), 'dmy');
    expect(dateOrder(['02/28/2026', '03/03/2026']), 'mdy');
    expect(dateOrder(['13/02/2026', '02/13/2026']), 'mixed');
    const unsure = 'Animal,Group,Date,Volume\nm1,V,02/03/2026,100\nm1,V,05/03/2026,150\n';
    const q = read(unsure);
    expect(waiting(q), true);
    expect(q.questions[0].options.map(([, l]) => l).join('|'), 'Day first: 02/03/2026 is 2 March 2026|' +
      'Month first: 02/03/2026 is February 3, 2026');
    expect(read(unsure, {dates: 'dmy'}).rows.map((x) => x.day).join(), '0,3');
    expect(read(unsure, {dates: 'mdy'}).rows.map((x) => x.day).join(), '0,89');
  });

  test('day 0 from a randomization date column, per animal, against the first measurement', async () => {
    const csv = 'Animal,Group,Date,Randomization date,Volume\n' +
      'm1,V,2026-03-02,2026-03-01,100\nm1,V,2026-03-05,2026-03-01,150\nm2,C,2026-03-04,2026-03-03,90\n';
    const r = read(csv);
    expect(r.rows.map((x) => x.day).join(), '1,4,1');
    expect(studySummary(r).includes('Day 0 is each animal\'s date in column D ("Randomization date").'), true,
      lines(r));
    const one = read(csv.replace('2026-03-03,90', '2026-03-01,90'));
    expect(one.rows.map((x) => x.day).join(), '1,4,3');
    expect(studySummary(one).includes('Day 0 is the date in column D ("Randomization date"), 2026-03-01.'), true);
    expect(read(csv.replace(/,2026-03-01,150/, ',2026-03-02,150')).error!.startsWith('Not read: rows 2 and 3 give ' +
      'one animal two different dates'), true);
  });

  test('two volume columns: the user picks one, never a guess', async () => {
    const csv = 'Animal,Group,Day,Volume,TV\nm1,V,0,100,101\nm2,C,0,90,91\n';
    const r = read(csv);
    expect(waiting(r), true);
    expect(r.questions[0].text, 'Which column is the tumor volume?');
    expect(r.questions[0].options.map(([, l]) => l).join(), 'Volume,TV');
    const picked = read(csv, {'role:volume': '4'});
    expect(picked.rows[0].volume, 101);
    expect(lines(picked).includes('Also kept: Volume (as given).'), true, lines(picked));
  });

  test('blinded: counts only, and no group name in any line or question', async () => {
    const csv = 'Animal,Group,Date,Length,Width,Comment\nm1,Alpha,2026-03-02,10,8,ok\nm2,Beta,2026-03-02,9,7,\n' +
      'Mean,,2026-03-02,9.5,7.5,\nm3,Gamma,,8,6,\n';
    const r = read(csv);
    const all = [...studySummary(r, true), ...r.questions.flatMap((q) => [q.text, ...q.options.map(([, l]) => l)])];
    expect(all.some((l) => /Alpha|Beta|Gamma/.test(l)), false, all.join(' | '));
    expect(studySummary(r, true)[0], 'Read 2 groups: 1 and 1 animals, days 0–0.');
    expect(r.ignored.join('|'), 'row 4, a summary row (Mean)|row 5, with no day');
  });

  test('every measurement cell placed, or left out with its reason: else nothing is read', async () => {
    const r = read('Animal,Group,Day,Volume\nm1,V,0,100\nm1,V,3,dead\n,V,7,50\nm2,C,0,90,7\n');
    expect(r.error, null);
    expect(r.ignored.join('|'), 'row 4, with no animal|1 cell that is not a plain number (D3)|column E, with no ' +
      'heading');
    expect(r.rows.length, 2);
    expect(read('Animal,Group,Day,Volume\nm1,V,0,100\nm1,V,0,110\n').error, 'Not read: rows 2 and 3 are one animal ' +
      'on one day; keep one row per animal and day.');
    expect(read('Animal,Group,Day,Volume\nm1,V,x,100\n').error!.startsWith('Not read: column C ("Day") holds ' +
      'neither days'), true);
    expect(read('Foo,Bar\n1,2\n').error!.startsWith('Not read: the columns are not named as a study\'s.'), true);
  });

  test('a wide table: read as before, its volumes\' unit assumed and said', async () => {
    const r = read('Animal,Group,Day 0,Day 3\nV1,Vehicle,37.9,51.8\nV2,,41.2,55.0\nC1,Cpd-A,39.5,45.0\n');
    expect(r.shape, 'by day');
    expect(lines(r), 'Read 2 groups: Vehicle 2 animals, Cpd-A 1 animal, days 0–3. | Volumes as given. | No unit in ' +
      'the column names: assumed mm³.');
    expect(lines(read('Animal,Group,Day 0 (mm3),Day 3 (mm3)\nV1,Vehicle,37.9,51.8\n')).includes('assumed'), false);
    expect(read('Animal,Group,Baseline,D3,D7\nm1,V,100,150,160\n').error!
      .startsWith('Not read: column C ("Baseline")'), true);
  });
  test('numbers under an empty heading are listed, wherever the column is', async () => {
    const r = read('Animal,Group,Day,Volume,,Note\n1,A,0,100,55,x\n2,A,0,120,66,y\n');
    expect(r.error, null);
    expect(r.ignored.join('|'), 'column E, with no heading');
    expect(studySummary(r).includes('Not read: column E, with no heading.'), true, lines(r));
  });

  test('an animal id in two groups: only the shared ids renamed, never by the group\'s name', async () => {
    const csv = 'Animal,Group,Day,Volume\n1,Alpha,0,100\n2,Alpha,0,110\n1,Beta,0,120\n3,Beta,0,90\n';
    const r = read(csv);
    expect(r.rows.map((x) => x.animal).join('|'), '1 (g1)|2|1 (g2)|3');
    expect(studySummary(r).includes('Animal ids repeat across groups: read as different animals (1 id).'), true);
    expect(r.rows.some((x) => /Alpha|Beta/.test(x.animal)), false);
    expect(studySummary(r, true).some((l) => /Alpha|Beta/.test(l)), false, studySummary(r, true).join('|'));
    expect(r.groups.map(([, n]) => n).join(), '2,2');
  });

  test('length and width but no volume, and no formula: left out and said', async () => {
    const r = read('Animal,Group,Day,Length,Width,Volume\n1,A,0,10,8,500\n1,A,3,12,9,\n');
    expect(r.ignored.join('|'), 'row 3, length and width but no volume, and no formula matches the volumes given');
    expect(r.rows[1].volume, null);
  });

  test('a day column that is not headed as a day says so; weeks are multiplied by 7; visit is no day', async () => {
    const t = read('Animal,Group,Time,Volume\nm1,V,0,100\nm1,V,3,150\n');
    expect(t.notes.includes('Days read from column C ("Time").'), true, lines(t));
    const d = read('Animal,Group,Day,Volume\nm1,V,0,100\nm1,V,3,150\n');
    expect(d.notes.some((n) => n.startsWith('Days read')), false);
    const w = read('Animal,Group,Week,Volume\nm1,V,0,100\nm1,V,2,150\n');
    expect(w.rows.map((x) => x.day).join(), '0,14');
    expect(w.notes.includes('Weeks in column C ("Week") read as days (× 7).'), true, lines(w));
    expect(waiting(read('Animal,Group,Visit,Volume\nm1,V,0,100\nm1,V,2,150\n')), true);
  });

  test('the height formulas first, and a formula chosen only when it beats the others', async () => {
    const sets = [[10, 8, 6], [12, 9, 7], [11, 7, 5], [11.5, 7.5, 6], [9, 6, 5]]
      .map(([l, w, h]) => ({l, w, h, v: FORMULAS['pi/6*LWH'].v(l, w, h)}));
    expect(deduceFormula(sets), 'pi/6*LWH');
    // the same volumes without a height column: the width-only formulas, none matching
    expect(deduceFormula(sets.map((x) => ({...x, h: null}))), null);
  });

  test('missing volumes filled from the deduced formula are counted', async () => {
    const csv = caliper('Animal,Group,Day,Length,Width,Volume', false, 'LW2/2').replace(/,323.44\n/, ',\n');
    const r = read(csv);
    expect(r.error, null);
    expect(studySummary(r)[1], 'Volumes match length × width² / 2; 1 missing volume computed with it.');
    expect(near(r.rows[3].volume, 11.5 * 7.5 * 7.5 / 2), true);
  });

  test('a day above 1000 is a date: the table is refused', async () => {
    expect(read('Animal,Group,Day,Volume\nm1,V,20260301,100\nm1,V,20260304,150\n').error, 'Not read: column C ' +
      '("Day") holds 20260301, not a study day: write dates as 2026-03-01.');
    expect(read('Animal,Group,Day 0,20260301\nV1,Vehicle,37.9,51.8\n').error!.startsWith('Not read: column D ' +
      '("20260301") holds 20260301, not a study day'), true);
  });

  test('a randomization date after the first measurement: negative days, said', async () => {
    const r = read('Animal,Group,Date,Randomization date,Volume\nm1,V,2026-03-01,2026-03-03,100\n' +
      'm1,V,2026-03-05,2026-03-03,150\n');
    expect(r.rows.map((x) => x.day).join(), '-2,2');
    expect(r.notes.includes('Measurements before day 0 are at negative days.'), true, lines(r));
  });

  test('a start-date column with no dates is not the start: kept as an extra column', async () => {
    const r = read('Animal,Group,Date,Randomized,Volume\nm1,V,2026-03-01,Yes,100\nm1,V,2026-03-05,Yes,150\n');
    expect(r.error, null);
    expect(r.rows.map((x) => x.day).join(), '0,4');
    expect(r.extra.map(([n]) => n).join(), 'Date,Randomized');
  });

  test('body weights that are very small for grams: the unit is asked', async () => {
    const csv = 'Animal,Group,Day,Volume (mm3),BW\nm1,V,0,100,0.02\nm2,C,0,120,0.021\n';
    const q = read(csv);
    expect(waiting(q), true);
    expect(q.questions[0].text, 'The body weights are very small for grams. In which unit are they?');
    expect(q.questions[0].options.map(([, l]) => l).join('|'), 'g, as they are|kg: multiplied by 1000');
    expect(read(csv, {'unit:weight': 'kg'}).rows[0].weight, 20);
    expect(read(csv, {'unit:weight': 'g'}).rows[0].weight, 0.02);
  });

  test('rows labelled as a summary with a group after the word are left out, said by the word alone', async () => {
    for (const label of ['Mean Vehicle', 'n=5 Vehicle', 'Mean, Vehicle']) {
      const r = read(`Animal,Group,Day,Volume\n"${label}",Vehicle,0,105\n1,Vehicle,0,100\n2,Vehicle,0,110\n`);
      expect(r.error, null, label);
      expect(r.rows.map((x) => x.animal).join(), '1,2', label);
      expect(r.ignored.length, 1, label);
      expect(/Vehicle/.test(r.ignored[0]), false, r.ignored[0]);
      expect(/^row 2, a summary row \((Mean|n)\)$/.test(r.ignored[0]), true, r.ignored[0]);
    }
    // a by-day table with the same label is refused, naming no group
    const w = read('Animal,Group,Day 0,Day 7\n1,Vehicle,100,200\n"Mean, Vehicle",,105,205\n');
    expect(w.error, 'Not read: row 3 looks like a summary row; leave its animal cell blank, or write only Mean (or ' +
      'n, SD).');
    // N-acetyl and SD-1 are names
    const names = read('Animal,Group,Day,Volume\nSD-1,N-acetyl,0,100\nSD-2,N-acetyl,0,110\n');
    expect(names.rows.length, 2);
    expect(names.ignored.length, 0);
  });

  test('a number in an extra column of a row that is left out is told, not dropped', async () => {
    const r = read('Animal,Group,Day,Volume,Notes\n1,A,0,100,\nMean,,0,105,n=5\n,,,,12\n');
    expect(r.error, null);
    expect(r.ignored.join(' | '), 'row 3, a summary row (Mean) | row 4, with no animal');
  });

  test('two columns that differ only in case are numbered once, in the reader, as the table will name them',
    async () => {
      const r = read('Animal,Group,Day,Volume,Cage,Cage,Comment,comment,Volume (as given)\n1,A,0,100,3,4,x,y,7\n' +
        '2,A,0,110,3,4,x,y,8\n');
      expect(r.extra.map(([n]) => n).join(), 'Cage,Cage (2),Comment,comment (2),Volume (as given)');
      expect(lines(r).includes('Also kept: Cage, Cage (2), Comment, comment (2), Volume (as given).'), true);
    });

  test('body weights with a median above 1000 are asked about: in mg?', async () => {
    const csv = 'Animal,Group,Day,Volume,Weight\n1,A,0,100,2100\n2,A,0,120,2200\n';
    const q = read(csv);
    expect(waiting(q), true);
    expect(q.questions[0].text, 'The body weights are very large for grams. In which unit are they?');
    expect(q.questions[0].options.map(([, l]) => l).join('|'), 'g, as they are|mg: divided by 1000');
    expect(read(csv, {'unit:weight': 'mg'}).rows[0].weight, 2.1);
    expect(read(csv, {'unit:weight': 'g'}).rows[0].weight, 2100);
  });

  test('an empty Day column beside a filled Date column: the dates are used', async () => {
    const r = read('Animal,Group,Date,Volume,Day\n1,A,2026-03-02,100,\n2,A,2026-03-02,110,\n' +
      '1,A,2026-03-09,150,\n');
    expect(r.error, null);
    expect(r.rows.map((x) => x.day).join(), '0,0,7');
    expect(r.notes.includes('Day 0 is the first measurement, 2026-03-02.'), true, lines(r));
  });

  test('cells that are not numbers: at most six named, the rest counted', async () => {
    const rows = Array.from({length: 14}, (_, i) => `m1,V,${i},x`).join('\n');
    const r = read(`Animal,Group,Day,Volume\n${rows}\nm2,V,0,50\n`);
    expect(r.ignored[0], '14 cells that are not a plain number (D2, D3, D4, D5, D6, D7, … and 8 more)');
  });
});

/** A sheet as PowerPack hands it over: every cell a string, a blank as null. */
const frame = (name: string, head: string[], rows: (string | null)[][]): SheetFrame => ({name, rowCount: rows.length,
  columns: {length: head.length, byIndex: (i) => ({name: head[i], getString: (r) => rows[r][i] as string})}});
const dated = (d: string) => `${d} 00:00:00 GMT+0100 (Central European Standard Time)`;

category('Excel sheets', () => {
  test('a sheet as cells: the names, then the rows, blanks as empty text', async () => {
    const cells = cellsOfFrame(frame('S', ['Animal', 'Group', 'Day'], [['m1', 'A', '0'], ['m2', null, '7']]));
    expect(JSON.stringify(cells), '[["Animal","Group","Day"],["m1","A","0"],["m2","","7"]]');
  });

  test('Excel dates read as days through the reader', async () => {
    const f = frame('S', ['Animal', 'Group', 'Date', 'Volume'], [['m1', 'A', dated('Mon Nov 02 2026'), '100'],
      ['m1', 'A', dated('Mon Nov 09 2026'), '150'], ['m2', 'B', dated('Mon Nov 02 2026'), '90'],
      ['m2', 'B', dated('Mon Nov 09 2026'), '95']]);
    const cells = cellsOfFrame(f);
    expect(cells[1][2], '2026-11-02');
    const r = readStudy(cells);
    expect(r.error === null, true);
    expect(JSON.stringify(r.rows.filter((x) => x.animal === 'm1').map((x) => x.day)), '[0,7]');
  });

  test('a sheet through CSV text reads the same, quoted cells included', async () => {
    const cells = [['Animal', 'Group', 'Day', 'Volume'], ['m1', 'A, high', '0', '100'], ['m1', 'A, high', '7', '150']];
    expect(JSON.stringify(cellsOf(csvOf(cells))), JSON.stringify(cells));
  });

  test('which sheets read as a study: none, one, several', async () => {
    const study = [['Animal', 'Group', 'Day', 'Volume'], ['m1', 'A', '0', '100'], ['m1', 'A', '7', '150']];
    const notes = [['Remark'], ['hello']];
    expect(studySheets([{name: 'Notes', cells: notes}]).length, 0);
    expect(studySheets([{name: 'Notes', cells: notes}, {name: 'Data', cells: study}]).join(), 'Data');
    expect(studySheets([{name: 'A', cells: study}, {name: 'B', cells: study}]).join(), 'A,B');
  });
});
