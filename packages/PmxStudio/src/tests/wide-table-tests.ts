import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {cellsOf, dayOf, readWide, wideSummary, type WideRead} from '../wide-table';
import {readFile} from '../open-study';

/** The pattern our own extractor broke on: each group's name over its first animal only (merged cells, blanks to
 *  the right), groups of 4, 3 and 5, and the last group's last animals (m11, m12) far right of its name. */
const BY_ANIMAL = ',Vehicle,,,,Cpd-1,,,Cpd-2,,,,\n' +
  'Day,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12\n' +
  '0,100,101,99,102,98,97,103,100,99,101,102,98\n' +
  '3,150,151,149,152,120,118,121,110,111,109,112,108\n' +
  '7,220,230,210,240,140,,139,115,117,116,118,119\n';

const read = (csv: string, answers = {}): WideRead => readWide(cellsOf(csv), answers)!;
const valuesIn = (csv: string) => cellsOf(csv).flat().filter((s) => /^\d+(\.\d+)?$/.test(s)).length;

category('Wide tables', () => {
  test('a column per animal, each group named over its first animal only: every animal read and counted', async () => {
    const w = read(BY_ANIMAL);
    expect(w.shape, 'by animal');
    expect(w.question, null);
    expect(w.error, null);
    expect(JSON.stringify(w.groups), JSON.stringify([['Vehicle', 4], ['Cpd-1', 3], ['Cpd-2', 5]]));
    expect(wideSummary(w)[0], 'Read 3 groups: Vehicle 4 animals, Cpd-1 3 animals, Cpd-2 5 animals.');
    const last = w.rows.filter((r) => r.animal === 'm11' || r.animal === 'm12');
    expect(last.length, 6);
    expect(last.every((r) => r.group === 'Cpd-2'), true);
    // a blank cell is a missing measurement, nothing else: 12 animals × 3 days, less m6 on day 7
    expect(w.rows.length, 35);
    expect(w.rows.find((r) => r.animal === 'm12' && r.day === 7)!.volume, 119);
  });

  test('the same table with trailing blanks cut from its header rows still reads the last animals', async () => {
    const cut = BY_ANIMAL.replace(',Cpd-2,,,,\n', ',Cpd-2\n').replace(',m11,m12\n', ',m11\n');
    const w = read(cut);
    expect(JSON.stringify(w.groups), JSON.stringify([['Vehicle', 4], ['Cpd-1', 3], ['Cpd-2', 5]]));
    expect(w.notes.some((n) => n.startsWith('No animal id over column M')), true);
  });

  test('a column per day, the group in a column: one row per animal and day', async () => {
    const csv = 'Animal,Group,D0,D3,Day 7\nm1,Vehicle,100,150,220\nm2,Vehicle,101,151,\nm3,Cpd,99,120,130\n';
    const w = read(csv);
    expect(w.shape, 'by day');
    expect(JSON.stringify(w.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
    expect(w.rows.length, valuesIn(csv));
    expect(w.rows.filter((r) => r.animal === 'm3').map((r) => r.day).join(), '0,3,7');
  });

  test('a column per day, the group in merged cells beside the animals or in a row over them', async () => {
    const merged = read(',,D0,D3\nVehicle,m1,100,150\n,m2,101,151\nCpd,m3,99,120\n,m4,98,110\n');
    expect(JSON.stringify(merged.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 2]]));
    const titled = read('Mouse,0,3,7\nVehicle\n1,100,150,220\n2,101,151,230\nCpd\n3,99,120,130\n');
    expect(JSON.stringify(titled.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
  });

  test('body weight is read from a block of its own, by its title', async () => {
    const w = read(',,Tumor volume,,Body weight,\nAnimal,Group,D0,D3,D0,D3\nm1,V,100,150,20,21\nm2,C,90,100,19,18\n');
    const m1 = w.rows.find((r) => r.animal === 'm1' && r.day === 3)!;
    expect(m1.volume, 150);
    expect(m1.weight, 21);
    // two blocks of days with no titles: which is the volume is asked, not guessed
    const q = read('Animal,Group,D0,D3,D0,D3\nm1,V,100,150,20,21\n').question!;
    expect(q.id, 'block');
    expect(q.choices.join('|'), 'columns C–D|columns E–F');
  });

  test('measurements before the first group name are asked about, never guessed or dropped', async () => {
    const csv = ',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n3,1,2,150,151,120,110\n';
    const w = read(csv);
    expect(w.rows.length, 0);
    expect(w.question!.id, 'before');
    expect(w.question!.text, 'Columns B, C hold measurements but come before the first group name. Which group ' +
      'are they in?');
    expect(JSON.stringify(read(csv, {before: 0}).groups), JSON.stringify([['Vehicle', 4], ['Cpd', 2]]));
    expect(JSON.stringify(read(csv, {before: 1}).groups),
      JSON.stringify([['No group name', 2], ['Vehicle', 2], ['Cpd', 2]]));
    const out = read(csv, {before: 2});
    expect(JSON.stringify(out.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 2]]));
    expect(wideSummary(out)[1], 'Not read: columns B, C, before the first group name (as you chose).');
  });

  test('a header that could be animals or groups, and days that could run either way, are asked about', async () => {
    expect(read('Day,A1,A2,A3\n0,100,101,99\n3,150,151,149\n').question!.id, 'header');
    expect(JSON.stringify(read('Day,A1,A2,A3\n0,100,101,99\n3,150,151,149\n', {header: 0}).groups),
      JSON.stringify([['One group', 3]]));
    // an empty corner over numbers going up: days down the first column, or animals across the days?
    expect(read(',D0,D3\n1,100,150\n2,110,160\n').question!.id, 'shape');
    // group names alone, merged: the animals are numbered within each group, and the user is told
    const unnamed = read('Day,Vehicle,,Cpd,,\n0,100,101,99,98,97\n3,150,151,149,140,130\n');
    expect(JSON.stringify(unnamed.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 3]]));
    expect(unnamed.rows.some((r) => r.animal === 'Cpd 3'), true);
  });

  test('what cannot be read is refused or named: numbers, a second table below, days', async () => {
    expect(read('Animal,Group,D0,D3\nm1,V,"1,234",150\n').error!.startsWith('Nothing read: 1 value could not'),
      true);
    const below = read('Animal,Group,D0,D3\nm1,V,100,150\nm2,V,101,151\n\nAnimal,Group,D0,D3\nm1,V,20,21\n');
    expect(below.left.join(), 'row 6, a second table below the first');
    expect(dayOf('D7'), 7);
    expect(dayOf('Day 14'), 14);
    expect(dayOf('0'), 0);
    expect(dayOf('2024-03-01'), null);
  });

  test('a long table is left to the long reading', async () => {
    expect(readWide(cellsOf('Mouse ID,Arm,Day,Volume (mm3)\nm1,V,0,100\nm1,V,3,150\n')), null);
    expect(readWide(cellsOf('Group,Day,Volume\nV,0,100\nV,3,150\n')), null);
    // a day column the long reading does not know by name: still not taken for a wide table
    expect(readWide(cellsOf('Animal,Group,Timepoint,Volume\nm1,V,0,100\nm1,V,3,150\nm2,V,0,90\n')), null);
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
    expect(wideSummary(r.wide!.read)[0], 'Read 3 groups: Vehicle 4 animals, Cpd-1 3 animals, Cpd-2 5 animals.');
    const asked = readFile('cro.csv', ',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n' +
      '3,1,2,150,151,120,110\n');
    expect(asked.kind === 'table' && asked.wide!.read.question!.id, 'before');
  });

  test('never silent: every number in a fixture is read, or sits in a column or row named as not read', async () => {
    let checked = 0;
    for (const [name, csv, given] of FIXTURES) {
      // a question answered with its first choice, as a user might: nothing goes missing whatever the answer
      const answers = {...given};
      let w = readWide(cellsOf(csv), answers);
      for (let k = 0; k < 4 && w?.question; k++)
        w = readWide(cellsOf(csv), Object.assign(answers, {[w.question.id]: 0}));
      if (!w || w.question || w.error)
        continue;
      checked++;
      expect(unaccounted(csv, w), 0, name);
    }
    expect(checked >= 25, true, `${checked} fixtures read`);
  });

  test('a day measured twice does not end the days: the row is named, every animal read', async () => {
    const w = read(',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n3,155,156,125\n7,160,170,130\n');
    expect(JSON.stringify(w.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
    expect(w.rows.filter((r) => r.animal === 'm1').map((r) => r.day).join(), '0,3,7');
    expect(w.left.join(), 'row 5, day 3 a second time');
  });

  test('body weight beside the volumes, a column per animal, is refused, not read as a group', async () => {
    const w = read(',Volume,,Weight,\nDay,m1,m2,m1,m2\n0,100,101,20,21\n3,150,151,19,20\n');
    expect(w.error, 'Nothing read: the names over the animals ("Volume", "Weight") are measurements, not groups. ' +
      'Open the tumor volumes on their own, each group\'s name over its first animal');
  });

  test('a row with an animal and no numbers is not taken for a group name', async () => {
    for (const second of ['2,,,', '2,dead,,']) {
      const w = read(`Mouse,0,3,7\nVehicle\n1,100,150,220\n${second}\n3,99,120,130\nCpd\n4,99,120,130\n`);
      expect(JSON.stringify(w.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]), second);
      expect(w.left.includes('row 4 (2), no numbers'), true, second);
    }
  });

  test('animal ids and group names written as numbers', async () => {
    const ids = read(',Vehicle,,Cpd\nDay,101,102,103\n0,100,101,99\n3,150,151,120\n7,160,170,130\n');
    expect(ids.shape, 'by animal');
    expect([...new Set(ids.rows.map((r) => r.animal))].join(), '101,102,103');
    expect(JSON.stringify(ids.groups), JSON.stringify([['Vehicle', 2], ['Cpd', 1]]));
    const groups = read(',1,,2\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n7,160,170,130\n');
    expect(JSON.stringify(groups.groups), JSON.stringify([['1', 2], ['2', 1]]));
  });

  test('an animal whose cells are all words is named as not read, in both layouts', async () => {
    const byAnimal = read(',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,dead,99\n3,150,dead,120\n');
    expect(byAnimal.left.includes('column C (m2), no numbers'), true);
    expect(byAnimal.left[0], 'cells that are not numbers: "dead" at C3, "dead" at C4');
    const byDay = read('Animal,Group,D0,D3\nm1,V,100,150\nm2,V,dead,dead\nm3,C,1,2\n');
    expect(byDay.left.includes('row 3 (m2), no numbers'), true);
  });

  test('a column that might hold the groups but is not headed so is asked about', async () => {
    const csv = 'Animal,Cage,D0,D3\nm1,C1,100,150\nm2,C1,101,151\nm3,C2,99,120\nm4,C2,98,110\n';
    const q = read(csv).question!;
    expect(q.id, 'group');
    expect(q.choices.join('|'), 'Column B ("Cage")|None: each group is named in a row over its animals');
    expect(JSON.stringify(read(csv, {group: 0}).groups), JSON.stringify([['C1', 2], ['C2', 2]]));
  });

  test('a summary of means is not taken for a wide table', async () => {
    expect(readWide(cellsOf('Group,Day,Mean,SEM\nV,0,100,5\nV,3,150,6\nV,7,200,9\n')), null);
    // group means across the top: asked, not read as animals unasked
    expect(read('Day,Vehicle,Cpd\n0,100,99\n3,150,120\n').question!.id, 'header');
  });

  test('the separator is the one the first lines use, outside quotes', async () => {
    expect(JSON.stringify(cellsOf('"a;b",c,d\n1;2,3,4\n')), JSON.stringify([['a;b', 'c', 'd'], ['1;2', '3', '4']]));
    expect(JSON.stringify(cellsOf('a;b;c\n1,5;2,5;3\n')), JSON.stringify([['a', 'b', 'c'], ['1,5', '2,5', '3']]));
  });

  test('blinded: the counts alone, in the table\'s order, and no group named in a question or a note', async () => {
    expect(wideSummary(read(BY_ANIMAL), true)[0], 'Read 3 groups: 4, 3 and 5 animals.');
    const q = read(',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n3,1,2,150,151,120,110\n')
      .question!;
    expect(q.blind!.choices[0], 'The group after them');
    const empty = read(',Vehicle,,Cpd,Spare\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n');
    expect(wideSummary(empty).includes('No measurements under "Spare".'), true);
    expect(wideSummary(empty, true).some((l) => l.includes('Spare') || l.includes('Vehicle')), false);
  });
});

/** The reviewer's probes: every layout read, asked about or refused; [name, csv, answers]. */
const FIXTURES: [string, string, {[id: string]: number}][] = [
  ['base', BY_ANIMAL, {}],
  ['ragged short rows', ',Vehicle,,Cpd\nDay,m1,m2,m3,m4\n0,100,101,99,98\n3,150,151\n7,220\n', {}],
  ['last group one animal', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', {}],
  ['group name in last col, no id', ',Vehicle,,Cpd\nDay,m1,m2,\n0,100,101,99\n3,150,151,120\n', {}],
  ['bom crlf', '\uFEFF,Vehicle,,Cpd\r\nDay,m1,m2,m3\r\n0,100,101,99\r\n3,150,151,120\r\n\r\n\r\n', {}],
  ['tabs', '\tVehicle\t\tCpd\nDay\tm1\tm2\tm3\n0\t100\t101\t99\n3\t150\t151\t120\n', {}],
  ['quoted comma group', ',"Vehicle, 10 mg",,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n', {}],
  ['NaN mid', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,NA,99\n3,150,dead,120\n7,160,170,abc\n', {}],
  ['day repeat', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n3,155,156,125\n7,160,170,130\n', {}],
  ['neg days', ',Vehicle,,Cpd\nDay,m1,m2,m3\n-1,100,101,99\n0,150,151,120\n3,160,170,130\n', {}],
  ['trailing blank cols', ',Vehicle,,Cpd,,,\nDay,m1,m2,m3,,,\n0,100,101,99,,,\n3,150,151,120,,,\n', {}],
  ['row first value at day', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,,\n3,150,151,120\n', {}],
  ['first data row single value', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,,\n3,,,120\n7,150,151,130\n', {}],
  ['days going back', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,1,2,3\n3,4,5,6\n7,7,8,9\n0,10,11,12\n', {}],
  ['before, own group', ',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n3,1,2,150,151,120,110\n',
    {before: 1}],
  ['before, left out', ',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n3,1,2,150,151,120,110\n',
    {before: 2}],
  ['by day ragged', 'Animal,Group,D0,D3,D7\nm1,V,100,150,220\nm2,V,101\nm3,C,99,120,130,140\n', {}],
  ['by day value beyond header', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,C,99,120,999\n', {}],
  ['by day NA', 'Animal,Group,D0,D3\nm1,V,NA,150\nm2,C,99,dead\n', {}],
  ['by day group row, blank animal', 'Mouse,0,3,7\nVehicle\n1,100,150,220\n2,,,\nCpd\n3,99,120,130\n', {}],
  ['by day animal with only text', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,dead,dead\nm3,C,1,2\n', {}],
  ['by day repeated day col', 'Animal,Group,D0,D3,D3,D7\nm1,V,100,150,151,200\nm2,C,99,120,121,140\n', {block: 0}],
  ['by day neg days', 'Animal,Group,D-1,D0,D3\nm1,V,100,150,151\nm2,C,99,120,121\n', {}],
  ['by day blank row inside', 'Animal,Group,D0,D3\nm1,V,100,150\n\nm2,C,99,120\n', {}],
  ['by day second table', 'Animal,Group,D0,D3\nm1,V,100,150\nm2,V,101,151\n\nAnimal,Group,D0,D3\nm1,V,20,21\n', {}],
  ['by day volume and weight', ',,Volume,,Weight,\nAnimal,Group,0,3,0,3\nm1,V,100,150,20,21\nm2,C,90,100,19,18\n', {}],
  ['numeric animal ids', ',Vehicle,,Cpd\nDay,101,102,103\n0,100,101,99\n3,150,151,120\n7,160,170,130\n', {}],
  ['numeric group names', ',1,,2\nDay,m1,m2,m3\n0,100,101,99\n3,150,151,120\n7,160,170,130\n', {}],
  ['group row, text-only animal',
    'Mouse,0,3,7\nVehicle\n1,100,150,220\n2,dead,,\n3,99,120,130\nCpd\n4,99,120,130\n', {}],
  ['by day merged group, text animal', 'Group,Animal,D0,D3\nV,m1,100,150\n,m2,101,151\nC,m3,99,120\n,m4,dead,\n', {}],
  ['by animal empty column between', ',Vehicle,,,Cpd\nDay,m1,m2,,m3\n0,100,101,,99\n3,150,151,,120\n', {}],
  ['by animal note left of days', 'Note,Day,m1,m2\nx,0,100,101\ny,3,150,151\n', {header: 0}],
  ['by animal all-text column', ',Vehicle,,Cpd\nDay,m1,m2,m3\n0,100,dead,99\n3,150,dead,120\n', {}],
  ['cage chosen', 'Animal,Cage,D0,D3\nm1,C1,100,150\nm2,C1,101,151\nm3,C2,99,120\nm4,C2,98,110\n', {group: 0}],
];

const COL = (s: string) => [...s].reduce((n, ch) => n * 26 + ch.charCodeAt(0) - 64, 0) - 1;
/** "C" or "C–F" (or "9–15") as the indexes it spans. */
const span = (x: string, of: (s: string) => number) => {
  const [a, b] = x.split('–').map(of);
  return Array.from({length: (b ?? a) - a + 1}, (_, k) => a + k);
};

/** The numbers of a fixture neither read nor in a column or row named as not read. Its own reading of the layout:
 *  by animal, the days down the column headed Day (else the first) under header rows; by day, the days across the
 *  first row that has two of them, labels to their left. */
function unaccounted(csv: string, w: WideRead): number {
  const cells = cellsOf(csv);
  const isNum = (s: string) => /^[+-]?\d[\d.,]*$/.test(s);
  const cols = new Set<number>();
  const rows = new Set<number>();
  for (const l of w.left) {
    const c = /^columns? ((?:[A-Z]+(?:–[A-Z]+)?(?:, )?)+)/.exec(l);
    c?.[1].split(', ').filter((x) => x).forEach((x) => span(x, COL)
      .forEach((k) => cols.add(k)));
    const r = /^rows? ((?:\d+(?:–\d+)?(?:, )?)+)/.exec(l);
    r?.[1].split(', ').filter((x) => x).forEach((x) => span(x, (s) => +s - 1)
      .forEach((k) => rows.add(k)));
  }
  let structural: (r: number, c: number) => boolean;
  if (w.shape === 'by animal') {
    const headed = cells.flatMap((row) => row.map((s, c) => /^day$/i.test(s) ? c : -1)).find((c) => c >= 0);
    const day = headed ?? Math.min(...cells.map((row) => row.findIndex((s) => s !== '')).filter((c) => c >= 0));
    const start = cells.findIndex((row) => dayOf(row[day] ?? '') !== null && row.some((s, c) => c > day && isNum(s)));
    structural = (r, c) => r < start || c === day;
  } else {
    const header = cells.findIndex((row) => row.filter((s, c) => c > 0 && dayOf(s) !== null).length >= 2);
    const first = cells[header].findIndex((s, c) => c > 0 && dayOf(s) !== null);
    structural = (r, c) => r <= header || c < first;
  }
  let numbers = 0;
  cells.forEach((row, r) => row.forEach((s, c) => {
    if (isNum(s) && !structural(r, c) && !rows.has(r) && !cols.has(c))
      numbers++;
  }));
  const read = w.rows.reduce((n, x) => n + (x.volume !== null ? 1 : 0) + (x.weight !== null ? 1 : 0), 0);
  return numbers - read;
}
