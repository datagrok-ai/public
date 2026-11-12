import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {SHAPES, cellsOf, dayOf, readWide, unplaced, wideSummary, type WideRead} from '../wide-table';
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

/** The numbers of a read table that are neither the source of exactly one row nor in a row or column left out with
 *  its reason, by position. The test's own reading of the layout: by animal, the rows down to the first one with a
 *  day in column A, and column A, are labels; by day, the first row and the columns left of the first day are. */
function unaccounted(csv: string, w: WideRead): string[] {
  const cells = cellsOf(csv);
  const isNum = (s: string) => /^[+-]?\d[\d.,]*$/.test(s);
  const rows = new Set(w.ignored.flatMap((l) => /^row (\d+)/.exec(l)?.[1] ?? []).map((r) => +r - 1));
  const cols = new Set(w.ignored.flatMap((l) => /^column ([A-Z]+)/.exec(l)?.[1] ?? [])
    .map((s) => [...s].reduce((n, ch) => n * 26 + ch.charCodeAt(0) - 64, 0) - 1));
  let label: (r: number, c: number) => boolean;
  if (w.shape === 'by animal') {
    const first = cells.findIndex((row) => dayOf(row[0] ?? '') !== null && row.slice(1).some(isNum));
    label = (r, c) => r < first || c === 0;
  } else {
    const top = cells.findIndex((row) => row.some((s) => s !== ''));
    const first = cells[top].findIndex((s) => s !== '' && dayOf(s) !== null);
    label = (r, c) => r <= top || c < first;
  }
  const from = new Map<string, number>();
  w.rows.forEach((x) => from.set(`${x.from[0]},${x.from[1]}`, (from.get(`${x.from[0]},${x.from[1]}`) ?? 0) + 1));
  const out: string[] = [];
  cells.forEach((row, r) => row.forEach((s, c) => {
    const n = from.get(`${r},${c}`) ?? 0;
    if (n > 1 || n === 1 && !isNum(s) || n === 0 && isNum(s) && !label(r, c) && !rows.has(r) && !cols.has(c))
      out.push(`${String.fromCharCode(65 + c)}${r + 1}=${s}`);
  }));
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
    const w = read(',Vehicle,,Cpd\nDay,101,102,103\nD0,100,101,99\nDay 3 (d),150,151,120\n');
    expect(groups(w), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
    expect([...new Set(w.rows.map((r) => r.animal))].join(), '101,102,103');
    expect(w.rows.map((r) => r.day).join(), '0,0,0,3,3,3');
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
    expect(byDay.ignored.join('|'), 'row 5 ("Mean"), a summary row|row 6 ("n"), a summary row');
    const tail = read('Animal,Group,D0,D3,\nm1,V,100,150,7\nm2,C,99,120,\n');
    expect(tail.ignored.join(), 'column E, with no heading, after the days');
    expect(wideSummary(tail)[1], 'Not read: column E, with no heading, after the days.');
    const byAnimal = read(BY_ANIMAL + 'SD,5,6,7,8,9,10,11,12,13,14,15,16\n');
    expect(byAnimal.ignored.join(), 'row 6 ("SD"), a summary row');
    expect(unaccounted(BY_ANIMAL + 'SD,5,6,7,8,9,10,11,12,13,14,15,16\n', byAnimal).join(), '');
  });

  test('anything else is refused with the two shapes that are read', async () => {
    expect(read(',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n3,155,156,125\n').error,
      `Not read: day 3 (row 5) does not come after the day above it. ${SHAPES}`);
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
      'Read 2 groups: 2 and 1 animals, days 0–7.|Not read: row 5 ("Mean"), a summary row.');
  });

  test('the separator is the one the first lines use; a comma is a decimal only between semicolons', async () => {
    expect(JSON.stringify(cellsOf('"a;b",c,d\n1;2,3,4\n')), JSON.stringify([['a;b', 'c', 'd'], ['1;2', '3', '4']]));
    const semi = 'Animal;Group;D0;D3\nm1;V;100,512;150,25\nm2;C;90,1;95,75\n';
    expect(readWide(cellsOf(semi), true)!.rows[0].volume, 100.512);
    expect(readWide(cellsOf(semi.replace(/;/g, '\t')), false)!.error!.startsWith('Not read: cell C2 ("100,512")'),
      true);
  });

  test('every reviewer input is read with each number placed, refused, or left to the long reading', async () => {
    for (const [name, csv, outcome] of FIXTURES) {
      const w = readWide(cellsOf(csv));
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
/** The reviewers' inputs (p1–p4) and what becomes of each: read (every number placed), refused, or left to the long
 *  reading as a table of one row per measurement. */
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
  ['6b cage + group', 'Animal,Group,Cage,D0,D3\nm1,V,C1,100,150\nm2,V,C1,101,151\nm3,T,C2,99,120\n', 'read'],
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
  ['N5b decimal comma 3 digits', 'Animal;Group;D0;D3\nm1;V;100,512;150,250\nm2;C;90,1;95,7\n', 'refused'],
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
];
/* eslint-enable max-len */
