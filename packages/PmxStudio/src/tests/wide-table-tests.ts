import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SHAPES, cellsOf, dayOf, readWide, separatorOf, unplaced, wideSummary, type WideRead} from '../wide-table';
import {readFile} from '../open-study';

/** The pattern our own extractor broke on: each group's name over its first animal only (merged cells, blanks to
 *  the right), groups of 4, 3 and 5, and the last group's last animals (m11, m12) far right of its name. */
const BY_ANIMAL = ',Vehicle,,,,Cpd-1,,,Cpd-2,,,,\n' +
  'Day,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12\n' +
  '0,100,101,99,102,98,97,103,100,99,101,102,98\n' +
  '3,150,151,149,152,120,118,121,110,111,109,112,108\n' +
  '7,220,230,210,240,140,,139,115,117,116,118,119\n';
const BY_DAY = 'Animal,Group,Day 0,Day 3,Day 7\nV1,Vehicle,37.9,51.8,88.5\nV2,,41.2,55.0,90.1\nC1,Cpd-A,39.5,45.0,\n';

const read = (csv: string): WideRead => readWide(cellsOf(csv))!;
const groups = (w: WideRead) => JSON.stringify(w.groups);

/** The cells with a digit in a read table that are neither the source of exactly one row nor in a row or column
 *  left out with its reason, and the rows whose day, animal or group is not the one at their cell's place. The
 *  test's own reading of the layout, none of the reader's. By animal, the first two filled rows (groups, ids) and
 *  column A are labels, the day is column A's number, the animal the id above and the group the nearest name to the
 *  left in the first row. By day, the first filled row and the columns left of its first cell with a digit are
 *  labels, and the day is that header's number. */
function unaccounted(csv: string, w: WideRead): string[] {
  const cells = cellsOf(csv);
  const filled = cells.map((_, r) => r).filter((r) => cells[r].some((s) => s !== ''));
  const rows = new Set(w.ignored.flatMap((l) => /^row (\d+)/.exec(l)?.[1] ?? []).map((r) => +r - 1));
  const cols = new Set(w.ignored.flatMap((l) => /^column ([A-Z]+)/.exec(l)?.[1] ?? [])
    .map((s) => [...s].reduce((n, ch) => n * 26 + ch.charCodeAt(0) - 64, 0) - 1));
  const numberIn = (s: string) => Number(/-?\d+(\.\d+)?/.exec(s)?.[0]);
  const byAnimal = w.shape === 'by animal';
  const [top, ids] = filled;
  const first = byAnimal ? 1 : cells[top].findIndex((s) => /\d/.test(s));
  const label = (r: number, c: number) => r <= (byAnimal ? ids : top) || c < first;
  const from = new Map<string, number>();
  w.rows.forEach((x) => from.set(`${x.from[0]},${x.from[1]}`, (from.get(`${x.from[0]},${x.from[1]}`) ?? 0) + 1));
  const out: string[] = [];
  cells.forEach((row, r) => row.forEach((s, c) => {
    const n = from.get(`${r},${c}`) ?? 0;
    if (n > 1 || n === 1 && !/\d/.test(s) || n === 0 && /\d/.test(s) && !label(r, c) && !rows.has(r) && !cols.has(c))
      out.push(`${String.fromCharCode(65 + c)}${r + 1}=${s}`);
  }));
  for (const x of w.rows) {
    const [r, c] = x.from;
    let group = '';
    for (let k = c; byAnimal && k >= 1 && !group; k--)
      group = cells[top][k] ?? '';
    if (x.day !== numberIn(byAnimal ? cells[r][0] : cells[top][c]) ||
      byAnimal && (x.group !== group || !x.animal.endsWith(cells[ids][c])))
      out.push(`${String.fromCharCode(65 + c)}${r + 1} placed as ${x.animal}, ${x.group}, day ${x.day}`);
  }
  return out;
}

category('Wide tables', () => {
  test('a column per animal, each group named over its first animal only: every animal read and counted', async () => {
    const w = read(BY_ANIMAL);
    expect(w.error, null);
    expect(w.shape, 'by animal');
    expect(groups(w), JSON.stringify([['Vehicle', 4], ['Cpd-1', 3], ['Cpd-2', 5]]));
    expect(wideSummary(w)[0], 'Read 3 groups: Vehicle 4 animals, Cpd-1 3 animals, Cpd-2 5 animals, days 0–7.');
    const last = w.rows.filter((r) => r.animal === 'm11' || r.animal === 'm12');
    expect(last.length, 6);
    expect(last.every((r) => r.group === 'Cpd-2'), true);
    // a blank cell is a missing measurement, nothing else: 12 animals × 3 days, less m6 on day 7
    expect(w.rows.length, 35);
    expect(w.rows.find((r) => r.animal === 'm12' && r.day === 7)!.volume, 119);
    expect(unaccounted(BY_ANIMAL, w).join(), '');
  });

  test('a column per animal with numeric ids, and days written D7 or with a unit', async () => {
    const w = read(',Vehicle,,Cpd\nDay,101,102,103\nD0,100,101,99\nDay 3 (mm3),150,151,120\n');
    expect(groups(w), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
    expect([...new Set(w.rows.map((r) => r.animal))].join(), '101,102,103');
    expect(w.rows.map((r) => r.day).join(), '0,0,0,3,3,3');
    // a weight's unit is not a day's; ids that are all numbers need "Day" before them; decimals are measurements
    expect(dayOf('Day 7 (g)'), null);
    expect(read(',Vehicle,,Cpd\n,101,102,103\n0,100,101,99\n').error!
      .startsWith('Not read: the animal ids are all numbers, so their row must start with "Day".'), true);
    expect(read(',Vehicle,,Cpd\n,37.9,41.2,39.5\n3,50,51,52\n').error!
      .startsWith('Not read: the animal id in column B is a number with decimals, a measurement.'), true);
  });

  test('a column per day, the group blank under its name: one row per animal and day', async () => {
    const w = read(BY_DAY);
    expect(w.shape, 'by day');
    expect(groups(w), JSON.stringify([['Vehicle', 2], ['Cpd-A', 1]]));
    expect(w.rows.length, 8);
    expect(w.rows.find((r) => r.animal === 'V2' && r.day === 3)!.group, 'Vehicle');
    expect(unaccounted(BY_DAY, w).join(), '');
  });

  test('rows of means and counts, and numbers in a column with no heading after the days, are listed', async () => {
    const byDay = read(BY_DAY + 'Mean,,39.5,50.6,89.3\nn,,3,3,2\n');
    expect(byDay.error, null);
    expect(byDay.ignored.join('|'), 'row 5, a summary row (Mean)|row 6, a summary row (n)');
    const tail = read('Animal,Group,D0,D3,\nm1,V,100,150,7\nm2,C,99,120,\n');
    expect(tail.ignored.join(), 'column E, with no heading, after the days');
    expect(wideSummary(tail)[1], 'Not read: column E, with no heading, after the days.');
    const byAnimal = read(BY_ANIMAL + 'SD,5,6,7,8,9,10,11,12,13,14,15,16\n');
    expect(byAnimal.ignored.join(), 'row 6, a summary row (SD)');
    // a summary is its word alone: N-acetyl is a group and SD-1 an animal, and the group is named again after one
    expect(groups(read('Animal,Group,D0,D3\nV1,Vehicle,1,2\nN1,N-acetyl,3,4\nN2,,5,6\n')),
      JSON.stringify([['Vehicle', 1], ['N-acetyl', 2]]));
    expect(groups(read('Animal,Group,D0,D3\nSD-1,Vehicle,1,2\nSD-2,,3,4\n')), JSON.stringify([['Vehicle', 2]]));
    expect(read('Animal,Group,D0,D3\nV1,Vehicle,1,2\nMean,,1,2\nV2,,3,4\n').error!
      .startsWith('Not read: row 4 has no group.'), true);
    expect(unaccounted(BY_ANIMAL + 'SD,5,6,7,8,9,10,11,12,13,14,15,16\n', byAnimal).join(), '');
  });

  test('anything else is refused with the two shapes that are read', async () => {
    expect(read(',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n3,155,156,125\n').error,
      `Not read: the day in row 5 does not come after the day above it. ${SHAPES}`);
    expect(read(',Vehicle,,Cpd,,Vehicle,\nDay,a,b,c,d,e,f\n0,1,2,3,4,5,6\n').error!
      .startsWith('Not read: the group named in column F starts again after another group.'), true);
    expect(SHAPES, 'PMX Studio reads two wide shapes. A column per animal: ",Vehicle,,Cpd-A" / "Day,V1,V2,C1" / ' +
      '"0,37.9,41.2,39.5" (each group\'s name over its first animal). A column per day: "Animal,Group,Day 0,Day 3" / ' +
      '"V1,Vehicle,37.9,51.8" / "V2,,41.2,55.0". Or open it as one row per measurement (animal, group, day, volume).');
    expect(read(',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n').error!
      .startsWith('Not read: the first animal (column B) has no group name above it.'), true);
    expect(read('Animal,Group,Arm,D0,D3\nm1,1,V,100,150\n').error!
      .startsWith('Not read: more than one column is headed as the groups ("Group", "Arm").'), true);
  });

  test('a number that could not be placed stops the read', async () => {
    const cells = cellsOf(BY_DAY);
    const w = read(BY_DAY);
    const label = (r: number, c: number) => r === 0 || c < 2;
    expect(unplaced(cells, w.rows, label, () => false), null);
    expect(unplaced(cells, w.rows.filter((x) => !(x.animal === 'V2' && x.day === 3)), label, () => false),
      'Nothing read: 1 number could not be placed (D3)');
  });

  test('blinded: the counts alone, in the table\'s order, and the days', async () => {
    expect(wideSummary(read(BY_ANIMAL), true)[0], 'Read 3 groups: 4, 3 and 5 animals, days 0–7.');
    expect(wideSummary(read(BY_DAY + 'Mean,,39.5,50.6,89.3\n'), true).join('|'),
      'Read 2 groups: 2 and 1 animals, days 0–7.|Not read: row 5, a summary row (Mean).');
  });

  test('blinded: no group\'s name in any line, read or refused', async () => {
    const inputs = [
      ',Alpha,,Beta,\nDay,a,b,c,d\n0,1,2,3,4\nMean Alpha,1.5,,3.5,\n',
      ',Alpha,,Beta,,Gamma\nDay,a,b,c,d\n0,1,2,3,4\n',
      ',Alpha,Beta,Alpha\nDay,a,b,c\n0,1,2,3\n',
      ',Alpha,,Body weight\nDay,a,b,c\n0,1,2,3\n',
      'Group,Animal,Day 0,Day 3\nAlpha,V1,1,2\n,V2,3,4\nAlpha,Mean,2,3\nBeta,,2,3\n',
      'Animal,Group,Day 0,Day 3\nV1,Alpha,1,2\nV2,,3,4\nMean Alpha,,2,3\n',
      'Animal,Group,Day 0,Day 3\nV1,Alpha,1,2\nV2,Beta,x Gamma,4\n',
      'Animal,Group,Day 0,Day 3\nV1,Alpha,1,2\nV2,Beta,3,4\n',
    ];
    for (const csv of inputs) {
      const w = read(csv);
      const lines = w.error ? [w.error.replace(SHAPES, '')] : wideSummary(w, true);
      expect(lines.some((l) => /Alpha|Beta|Gamma/.test(l)), false, lines.join(' | '));
    }
  });

  test('the separator is the one the first lines use; a comma is a decimal only between semicolons', async () => {
    expect(JSON.stringify(cellsOf('"a;b",c,d\n1;2,3,4\n')), JSON.stringify([['a;b', 'c', 'd'], ['1;2', '3', '4']]));
    const semi = 'Animal;Group;D0;D3\nm1;V;100,512;150,25\nm2;C;90,1;95,75\n';
    expect(readWide(cellsOf(semi), true)!.rows[0].volume, 100.512);
    expect(readWide(cellsOf(semi.replace(/;/g, '\t')), false)!.error!
      .startsWith('Not read: cell C2 is not a plain number'), true);
    // between semicolons, a dot may mark the thousands: 1.234 is refused
    expect(readWide(cellsOf('Animal;Group;D0;D3\nm1;V;1.234;2\n'), true)!.error!
      .startsWith('Not read: cell C2 is not a plain number'), true);
  });

  test('every reviewer input is read with each number placed, refused, or left to the long reading', async () => {
    for (const [name, csv, outcome] of FIXTURES) {
      const w = readWide(cellsOf(csv), separatorOf(csv) === ';');
      expect(w === null ? 'long' : w.error ? 'refused' : 'read', outcome, name);
      if (w && !w.error)
        expect(unaccounted(csv, w).join(), '', name);
    }
  });

  test('Open a study: a wide file opens as the long table, with its animals per group said', async () => {
    const r = readFile('cro.csv', BY_ANIMAL);
    expect(r.kind, 'table');
    if (r.kind !== 'table')
      return;
    expect(r.df.rowCount, 35);
    expect(r.df.columns.names().join(), 'animal,group,day,volume');
    expect(r.animals, 12);
    expect(r.groups, 3);
    expect(r.lastDay, 7);
    expect(wideSummary(r.wide!)[0], 'Read 3 groups: Vehicle 4 animals, Cpd-1 3 animals, Cpd-2 5 animals, days 0–7.');
    const refused = readFile('cro.csv', 'Animal,Group,Baseline,D3\nm1,V,100,150\n');
    expect(refused.kind === 'table' && refused.wide!.error!.startsWith('Not read: column C ("Baseline")'), true);
  });
});

/* eslint-disable max-len */
/** The reviewers' inputs (p1–p4; R3, the third review's) and what becomes of each: read (every number placed),
 *  refused, or left to the long reading as a table of one row per measurement. */
const FIXTURES: [string, string, 'read' | 'refused' | 'long'][] = [
  ['1 repeated day', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n3,155,156,125\n7,160,170,130\n10,170,180,140\n', 'refused'],
  ['2 weight block beside by animal', ',Volume,,Weight,\nDay,m1,m2,m1,m2\n0,100,101,20,21\n3,150,151,19,20\n', 'refused'],
  ['2b weight block group-named', ',Vehicle,,Vehicle BW,\nDay,m1,m2,m1,m2\n0,100,101,20,21\n3,150,151,19,20\n', 'refused'],
  ['3 group-row blank animal', 'Mouse,0,3,7\nVehicle\n1,100,150,220\n2,,,\n3,99,120,130\nCpd\n4,99,120,130\n', 'refused'],
  ['3b group-row word-only animal', 'Mouse,0,3,7\nVehicle\n1,100,150,220\n2,dead,,\n3,99,120,130\nCpd\n4,99,120,130\n', 'refused'],
  ['4 numeric ids', ',Vehicle,,Cpd\nDay,101,102,103\n0,100,101,99\n3,150,151,120\n7,160,170,130\n', 'read'],
  ['4b numeric group names', ',1,,2\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n7,160,170,130\n', 'read'],
  ['4c numeric group names by day', 'Animal,Group,D0,D3\n1,1,100,150\n2,1,101,151\n3,2,99,120\n', 'read'],
  ['5 text-only animal by animal', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,dead,99\n3,150,dead,120\n', 'refused'],
  ['5b text-only animal by day', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,dead,dead\nm3,C,1,2\n', 'refused'],
  ['6 cage column', 'Animal,Cage,D0,D3\nm1,C1,100,150\nm2,C1,101,151\nm3,C2,99,120\nm4,C2,98,110\n', 'refused'],
  ['6b cage + group', 'Animal,Group,Cage,D0,D3\nm1,V,C1,100,150\nm2,V,C1,101,151\nm3,T,C2,99,120\n', 'refused'],
  ['7 summary', 'Group,Day 0,Day 3,Day 7\nVehicle,100,150,200\nCpd,99,120,130\n', 'refused'],
  ['7b summary mean rows', ',Vehicle,Cpd\nDay,Mean,Mean\n0,100,99\n3,150,120\n', 'refused'],
  ['7c summary by animal-like', 'Day,Vehicle,Cpd\n0,100,99\n3,150,120\n7,200,130\n', 'refused'],
  ['8 multi-line quoted', ',"Vehicle\n10 mg/kg",,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', 'read'],
  ['8b separator: first lines text', 'Study X; report\nNotes; a; b\n,Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n7,1,2,3\n', 'long'],
  ['N1 two header rows of ids', ',Vehicle,,Cpd\nDay,m1,m2,m3\n,A,B,C\n0,100,101,99\n3,150,151,120\n', 'refused'],
  ['N2 group repeated over every animal', ',Vehicle,Vehicle,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', 'read'],
  ['N3 group header spans into empty separator col', ',Vehicle,,,Cpd\nDay,m1,m2,,m3\n0,100,101,,99\n3,150,151,,120\n', 'refused'],
  ['N3b separator col holds a stray number', ',Vehicle,,,Cpd\nDay,m1,m2,,m3\n0,100,101,5,99\n3,150,151,,120\n', 'refused'],
  ['N4 day cols with units', 'Animal,Group,Day 7 (mm3),Day 14 (mm3)\nm1,V,100,150\nm2,C,90,95\n', 'read'],
  ['N5 decimal comma semicolon', 'Animal;Group;D0;D3\nm1;V;100,5;150,25\nm2;C;90,1;95,75\n', 'read'],
  ['N5b decimal comma 3 digits', 'Animal;Group;D0;D3\nm1;V;100,512;150,250\nm2;C;90,1;95,7\n', 'read'],
  ['N6 days as D-1', 'Animal,Group,D-1,D0,D3\nm1,V,100,150,151\nm2,C,99,120,121\n', 'read'],
  ['N7 descending days by day', 'Animal,Group,D7,D3,D0\nm1,V,200,150,100\nm2,C,99,120,121\n', 'refused'],
  ['N7b descending days by animal', ',Vehicle,,Cpd\nDay,m1,m2,m3\n7,100,101,99\n3,150,151,120\n0,1,2,3\n', 'refused'],
  ['N8 totals row by day', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,101,151\nm3,C,99,120\nMean,,100,140.3\n', 'read'],
  ['N8b mean row by animal (bottom)', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\nMean,125,126,110\n', 'read'],
  ['N8c mean column by animal', ',Vehicle,,,Cpd\nDay,m1,m2,Mean,m3\n0,100,101,100.5,99\n3,150,151,150.5,120\n', 'refused'],
  ['N8d mean row inside group by day group rows', 'Mouse,0,3\nVehicle\n1,100,150\n2,101,151\nMean,100.5,150.5\nCpd\n3,99,120\n', 'refused'],
  ['N9 n= row', ',Vehicle,,Cpd\nn=,2,,1\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', 'refused'],
  ['N9b n= row by day', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,C,99,120\nn,,2,2\n', 'read'],
  ['N10 footnotes', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,C,99,120\n\n* animal 3 found dead day 5\nSource: CRO 2024\n', 'refused'],
  ['N10b footnote with number by animal', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n\nNote: 1 animal excluded,,,\n', 'refused'],
  ['N11 trailing commas excel', 'Animal,Group,D0,D3,,,\nm1,V,100,150,,,\nm2,C,99,120,,,\n,,,,,,\n,,,,,,\n', 'read'],
  ['N12 one group by day', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,99,120\n', 'read'],
  ['N12b one group by animal', ',Vehicle,\nDay,m1,m2\n0,100,101\n3,150,151\n', 'read'],
  ['N13 one animal by day', 'Animal,Group,D0,D3,D7\nm1,V,100,150,200\n', 'read'],
  ['N13b one animal by animal', ',Vehicle\nDay,m1\n0,100\n3,150\n', 'read'],
  ['N14 group name equals animal id', 'Animal,Group,D0,D3\nV,V,100,150\nm2,C,99,120\n', 'read'],
  ['N14b by animal group == id', ',m1,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', 'read'],
  ['N15 all NA column', 'Animal,Group,D0,D3,D7\nm1,V,100,150,NA\nm2,C,99,120,NA\n', 'read'],
  ['N16 blank lines between groups by day', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,101,151\n\n\nm3,C,99,120\n', 'read'],
  ['N16b blank lines between group-title blocks', 'Mouse,D0,D3\n\nVehicle\n1,100,150\n2,101,151\n\nCpd\n3,99,120\n', 'refused'],
  ['N17 dup ids across groups', 'Animal,Group,D0,D3\n1,Vehicle,100,150\n2,Vehicle,101,151\n1,Cpd,99,120\n2,Cpd,98,110\n', 'read'],
  ['N18 Treatment+Dose cols', 'Animal,Treatment,Dose,D0,D3\nm1,V,0,100,150\nm2,C,10,99,120\n', 'refused'],
  ['N19 day header in 2nd row, title row', 'Tumor volume (mm3),,,\nAnimal,Group,D0,D3\nm1,V,100,150\nm2,C,99,120\n', 'refused'],
  ['N20 days with "Day" prefix by animal col', ',Vehicle,,Cpd\nStudy day,m1,m2,m3\nDay 0,100,101,99\nDay 3,150,151,120\n', 'read'],
  ['N21 dates in first col', ',Vehicle,,Cpd\nDate,m1,m2,m3\n2024-01-01,100,101,99\n2024-01-04,150,151,120\n', 'long'],
  ['N22 negative / zero vol, "<LOD"', 'Animal,Group,D0,D3\nm1,V,0,<20\nm2,C,99,120\n', 'refused'],
  ['N23 group col then animal col numbers only', 'Group,Animal,D0,D3\n1,1,100,150\n1,2,101,151\n2,3,99,120\n', 'read'],
  ['N24 day header row repeated per group (stacked tables)', 'Vehicle,,\nAnimal,D0,D3\nm1,100,150\nm2,101,151\nCpd,,\nAnimal,D0,D3\nm3,99,120\n', 'refused'],
  ['N25 two group columns: Group and Arm', 'Animal,Group,Arm,D0,D3\nm1,1,V,100,150\nm2,2,C,99,120\n', 'refused'],
  ['N26 by animal value in day col left', 'Grp,Day,m1,m2\n,0,100,101\n,3,150,151\n', 'long'],
  ['N27 thousands', 'Animal,Group,D0,D3\nm1,V,"1,234.5",150\nm2,C,99,"1,200"\n', 'refused'],
  ['N28 sex col', 'Animal,Sex,Group,D0,D3\nm1,F,V,100,150\nm2,M,C,99,120\n', 'read'],
  ['N29 Time column header numbers only by animal, no Day word', ',Vehicle,,Cpd\n,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', 'read'],
  ['N29b numeric ids under a blank, no Day word', ',Vehicle,,Cpd\n,101,102,103\n0,100,101,99\n3,150,151,120\n', 'refused'],
  ['N30 by day plain number days with label "Mouse"', 'Mouse,Group,0,3,7\n1,V,100,150,200\n2,C,99,120,130\n', 'read'],
  ['N31 by day, unlabelled numeric header', ',,0,3,7\n1,V,100,150,200\n2,C,99,120,130\n', 'refused'],
  ['M1 text-id animal with no values under group titles', 'Mouse,0,3\nVehicle\nm1,100,150\nm2\nm3,99,120\nCpd\nm4,98,110\n', 'refused'],
  ['M1b same with trailing commas', 'Mouse,0,3\nVehicle,,\nm1,100,150\nm2,,\nm3,99,120\nCpd,,\nm4,98,110\n', 'refused'],
  ['M2 baseline col left of days', 'Animal,Group,Baseline,D3,D7\nm1,V,100,150,200\nm2,C,99,120,130\n', 'refused'],
  ['M2b BW single column left', 'Animal,Group,BW,D0,D3\nm1,V,20,100,150\nm2,C,21,99,120\n', 'refused'],
  ['M2c randomization volume col', 'Animal,Group,Rand vol,D0,D3\nm1,V,95,100,150\nm2,C,94,99,120\n', 'refused'],
  ['M3 cage col after group (unheaded group?)', 'Animal,Group,Cage,D0,D3\nm1,V,1,100,150\nm2,C,2,99,120\n', 'refused'],
  ['M4 group col blank first rows', 'Animal,Group,D0,D3\nm1,,100,150\nm2,V,101,151\nm3,C,99,120\n', 'refused'],
  ['M5 two group-ish cols unheaded', 'Animal,,,D0,D3\nm1,V,A,100,150\nm2,V,A,101,151\nm3,C,B,99,120\n', 'refused'],
  ['M6 day row with trailing text col (comments)', 'Animal,Group,D0,D3,Comment\nm1,V,100,150,ok\nm2,C,99,120,dead 5\n', 'refused'],
  ['M6b numeric comment col after days', 'Animal,Group,D0,D3,TGI%\nm1,V,100,150,\nm2,C,99,120,45\n', 'refused'],
  ['M7 by animal TGI row at bottom with day number', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n%TGI,,,40\n', 'refused'],
  ['M8 by animal: comment row between', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\nDosing,x,x,x\n3,150,151,120\n', 'refused'],
  ['L1 long basic', 'Animal,Group,Day,Volume\nm1,V,0,100\nm1,V,3,150\nm2,C,0,99\nm2,C,3,120\n', 'long'],
  ['L2 long +cage sex bw', 'Animal,Cage,Sex,Group,Day,Volume,Body weight\nm1,1,F,V,0,100,20\nm1,1,F,V,3,150,21\nm2,2,M,C,0,99,22\nm2,2,M,C,3,120,23\n', 'long'],
  ['L3 long day first', 'Day,Animal,Group,Volume\n0,m1,V,100\n0,m2,C,99\n3,m1,V,150\n3,m2,C,120\n', 'long'],
  ['L4 long day first numeric animals', 'Day,Mouse,Group,TV\n0,1,1,100\n0,2,2,99\n3,1,1,150\n3,2,2,120\n7,1,1,200\n7,2,2,130\n', 'long'],
  ['L5 long no Day word (Time)', 'ID,Arm,Time,Vol\n1,V,0,100\n1,V,3,150\n2,C,0,99\n2,C,3,120\n', 'long'],
  ['L6 long headers odd: Tag,Tx,DPI,TV', 'Tag,Tx,DPI,TV\nt1,V,0,100\nt1,V,3,150\nt2,C,0,99\nt2,C,3,120\n', 'long'],
  ['L7 long sorted by day, numeric everything', 'Day,Animal,Cage,Volume,Weight\n0,1,1,100,20\n0,2,1,101,21\n3,1,1,150,21\n3,2,1,151,22\n7,1,1,200,22\n7,2,1,201,23\n', 'long'],
  ['L8 long semicolon decimal comma', 'Animal;Group;Day;Volume\nm1;V;0;100,5\nm1;V;3;150,2\n', 'long'],
  ['L9 long one animal per day sorted by day lowercase', 'day,animal,group,volume\n0,1,V,100\n3,1,V,150\n7,1,V,200\n', 'long'],
  ['L10 long with Study day header + Date', 'Subject,Group,Date,Study Day,Tumor Volume (mm3)\nm1,V,2024-01-01,0,100\nm1,V,2024-01-04,3,150\n', 'long'],
  ['L11 long DPI no names first col days', 'DPI,Mouse,Group,TV\n0,m1,V,100\n3,m1,V,150\n7,m1,V,200\n0,m2,C,90\n3,m2,C,95\n', 'long'],
  ['L12 long day col named "Day" second, cage before', 'Cage,Day,Animal,Group,Volume\n1,0,m1,V,100\n1,3,m1,V,150\n1,7,m1,V,200\n2,0,m2,C,90\n', 'long'],
  ['L13 long, Days post implant + day ascending all', 'Days,Volume,Weight,Animal\n0,100,20,m1\n3,150,21,m1\n7,200,22,m1\n10,250,22,m1\n', 'long'],
  ['side-by-side group blocks by day', ',,Vehicle,,,Cpd-X,,\nAnimal,,D0,D3,D7,D0,D3,D7\nm1,,100,150,200,,,\nm2,,101,151,201,,,\nm3,,,,,99,120,130\n', 'refused'],
  ['N31 answer 1', ',,0,3,7\n1,V,100,150,200\n2,C,99,120,130\n', 'refused'],
  ['R3 S1 base', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4\n3,5,6,7,8\n', 'read'],
  ['R3 S1 noncontig repeat', ',Vehicle,Cpd,Vehicle\nDay,a,b,c\n0,1,2,3\n3,4,5,6\n', 'refused'],
  ['R3 S1 noncontig merged', ',Vehicle,,Cpd,,Vehicle,\nDay,a,b,c,d,e,f\n0,1,2,3,4,5,6\n', 'refused'],
  ['R3 S1 group over blank id', ',Veh,,Cpd\nDay,a,b,\n0,1,2,3\n', 'refused'],
  ['R3 S1 group over blank id mid', ',Veh,Cpd,\nDay,a,,c\n0,1,2,3\n', 'refused'],
  ['R3 S1 id row is data (col A blank)', ',Veh,,Cpd,\n,37.9,41.2,39.5,40.1\n3,50,51,52,53\n7,60,61,62,63\n', 'refused'],
  ['R3 S1 id row is data integers', ',Veh,,Cpd,\n,100,101,99,102\n3,150,151,149,152\n', 'refused'],
  ['R3 S1 summary mid', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4\nMean,1,2,3,4\n3,5,6,7,8\n', 'refused'],
  ['R3 S1 summary end blinded label', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4\nMean Veh,1.5,,3.5,\n', 'refused'],
  ['R3 S1 trailing seps', ',Veh,,Cpd,,\nDay,a,b,c,d,\n0,1,2,3,4,\n3,5,6,7,8,,\n', 'read'],
  ['R3 S1 ragged short', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2\n3,5,6,7,8\n', 'read'],
  ['R3 S1 ragged long', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4,99\n', 'refused'],
  ['R3 S1 dup day', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4\n0,5,6,7,8\n', 'refused'],
  ['R3 S1 D7 vs 7', ',Veh,,Cpd,\nDay,a,b,c,d\nD7,1,2,3,4\n7.0,5,6,7,8\n', 'refused'],
  ['R3 S1 BOM', '\uFEFF,Veh,,Cpd,\nDay,a,b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S1 CRLF', ',Veh,,Cpd,\r\nDay,a,b,c,d\r\n0,1,2,3,4\r\n', 'read'],
  ['R3 S1 quoted comma value', ',Veh,,Cpd,\nDay,a,b,c,d\n0,"37,9",2,3,4\n', 'read'],
  ['R3 S1 quoted 1,234', ',Veh,,Cpd,\nDay,a,b,c,d\n0,"1,234",2,3,4\n', 'refused'],
  ['R3 S1 quoted newline id', ',Veh,,Cpd,\nDay,"a\nx",b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S1 semicolon 1.234 thousands', ';Veh;;Cpd;\nDay;a;b;c;d\n0;1.234;2,5;3;4\n', 'refused'],
  ['R3 S1 semicolon 100,512', ';Veh;;Cpd;\nDay;a;b;c;d\n0;100,512;2,5;3;4\n', 'read'],
  ['R3 S1 group Volasertib', ',Vehicle,,Volasertib,\nDay,a,b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S1 beyond group blinded', ',Veh,,Cpd,,Secret\nDay,a,b,c,d\n0,1,2,3,4\n', 'refused'],
  ['R3 S1 dup ids across groups', ',Veh,,Cpd,\nDay,1,2,1,2\n0,1,2,3,4\n', 'read'],
  ['R3 S1 id named N-1', ',Veh,,Cpd,\nDay,N-1,b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S1 numbers 1e3', ',Veh,,Cpd,\nDay,a,b,c,d\n0,1e3,2,3,4\n', 'refused'],
  ['R3 S1 value with star', ',Veh,,Cpd,\nDay,a,b,c,d\n0,37.9*,2,3,4\n', 'refused'],
  ['R3 S1 day with (g) unit', ',Veh,,Cpd,\nDay,a,b,c,d\nDay 0 (mm3),1,2,3,4\nDay 7 (g),20,21,22,23\n', 'refused'],
  ['R3 S1 group repeated each animal', ',Veh,Veh,Cpd,Cpd\nDay,a,b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S1 group row has number group', ',1,,2,\nDay,a,b,c,d\n0,1,2,3,4\n', 'read'],
  ['R3 S2 base', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2\nV2,,3,4\nC1,Cpd,5,6\n', 'read'],
  ['R3 S2 N-acetyl carry', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2\nN1,N-acetyl,3,4\nN2,,5,6\n', 'read'],
  ['R3 S2 SD-1 id', 'Animal,Group,Day 0,Day 3\nSD-1,Vehicle,1,2\nSD-2,,3,4\n', 'read'],
  ['R3 S2 group-first summary blinded', 'Group,Animal,Day 0,Day 3\nVehicle,V1,1,2\n,V2,3,4\nVehicle,Mean,2,3\n', 'read'],
  ['R3 S2 mean label w group blinded', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2\nV2,,3,4\nMean Vehicle,,2,3\n', 'refused'],
  ['R3 S2 dup day header', 'Animal,Group,Day 7,D7\nV1,Vehicle,1,2\n', 'refused'],
  ['R3 S2 merged day label', 'Animal,Group,Day,,\n,,0,3,7\nV1,Vehicle,1,2,3\n', 'refused'],
  ['R3 S2 day col actually volumes', 'Animal,Group,100,200\nV1,Vehicle,1,2\n', 'read'],
  ['R3 S2 weight col with day header', 'Animal,Group,Day 0 (mm3),Day 7 (mm3),Day 14 (g)\nV1,Vehicle,1,2,20\n', 'refused'],
  ['R3 S2 notes col with number star', 'Animal,Group,Notes,Day 0,Day 3\nV1,Vehicle,37.9*,1,2\n', 'refused'],
  ['R3 S2 notes col 1e3', 'Animal,Group,Notes,Day 0,Day 3\nV1,Vehicle,1e3,1,2\n', 'refused'],
  ['R3 S2 blank hdr after days', 'Animal,Group,Day 0,Day 3,\nV1,Vehicle,1,2,99\n', 'read'],
  ['R3 S2 blank hdr between days', 'Animal,Group,Day 0,,Day 3\nV1,Vehicle,1,9,2\n', 'refused'],
  ['R3 S2 extra header ragged long row', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2,3,4\n', 'read'],
  ['R3 S2 title above', 'Study X\nAnimal,Group,Day 0,Day 3\nV1,Vehicle,1,2\n', 'refused'],
  ['R3 S2 two group cols', 'Animal,Group,Treatment,Day 0,Day 3\nV1,Vehicle,X,1,2\n', 'refused'],
  ['R3 S2 dose col numbers', 'Animal,Group,Dose,Day 0,Day 3\nV1,Vehicle,10,1,2\n', 'refused'],
  ['R3 S2 first row blank group', 'Animal,Group,Day 0,Day 3\nV1,,1,2\n', 'refused'],
  ['R3 S2 numeric day hdr + id', 'ID,Group,0,3,7\nV1,Vehicle,1,2,3\n', 'read'],
  ['R3 S2 summary mid then carry', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2\nMean,,1,2\nV2,,3,4\n', 'refused'],
  ['R3 S2 N group summary row', 'Animal,Group,Day 0,Day 3\nV1,Vehicle,1,2\nn,,1,1\n', 'read'],
  ['R3 S2 semicolon', 'Animal;Group;Day 0;Day 3\nV1;Vehicle;1,5;2\n', 'read'],
  ['R3 S2 semicolon 1.234', 'Animal;Group;Day 0;Day 3\nV1;Vehicle;1.234;2\n', 'refused'],
  ['R3 S2 tab', 'Animal\tGroup\tDay 0\tDay 3\nV1\tVehicle\t1\t2\n', 'read'],
  ['R3 S2 quoted sep in group', 'Animal,Group,Day 0,Day 3\nV1,"Cpd, 10 mg/kg",1,2\n', 'read'],
  ['R3 LONG base', 'Animal,Group,Day,Volume\nM1,Veh,0,100\nM1,Veh,3,120\n', 'long'],
  ['R3 LONG timepoint', 'Mouse,Arm,Timepoint,TV\nM1,Veh,0,100\nM1,Veh,3,120\n', 'long'],
  ['R3 LONG pandas index Week', ',Animal,Group,Week,Volume\n0,M1,Veh,0,100\n1,M1,Veh,1,120\n2,M2,Veh,0,90\n', 'long'],
  ['R3 LONG D-coded day', 'Mouse,Arm,Timepoint,TV\nM1,Veh,D0,100.5\nM1,Veh,D3,120\n', 'long'],
  ['R3 LONG visit', 'Subject,Treatment,Visit,Volume\nS1,Veh,Day 0,100\nS1,Veh,Day 3,120\n', 'long'],
  ['R3 LONG blank first header', ',Group,Day,Volume\nM1,Veh,0,100\n', 'long'],
  ['R3 LONG no header', 'M1,Veh,0,100\nM1,Veh,3,120\nM2,Veh,0,90\n', 'long'],
  ['R3 empty', '', 'long'],
  ['R3 one cell', '5', 'long'],
  ['R3 header only S2', 'Animal,Group,Day 0,Day 3\n', 'long'],
  ['R3 two rows only S1 hdr', ',Veh,,Cpd,\nDay,a,b,c,d\n', 'long'],
];
/* eslint-enable max-len */
