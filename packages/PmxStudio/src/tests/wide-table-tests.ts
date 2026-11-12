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
    // "Cage" over numbers going up: days down the first column, or cages across the days?
    expect(read('Cage,D0,D3\n1,100,150\n2,110,160\n').question!.id, 'shape');
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
    expect(r.wide!.summary[0], 'Read 3 groups: Vehicle 4 animals, Cpd-1 3 animals, Cpd-2 5 animals.');
    const asked = readFile('cro.csv', ',,,Vehicle,,Cpd,\nDay,m0,m00,m1,m2,m3,m4\n0,1,2,100,101,99,98\n' +
      '3,1,2,150,151,120,110\n');
    expect(asked.kind === 'table' && asked.wide!.question!.id, 'before');
  });
});
