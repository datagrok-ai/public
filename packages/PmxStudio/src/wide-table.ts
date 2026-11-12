/** The columns Analyze needs, by the names a study table usually gives them. */
export const ROLES: [string, RegExp][] = [
  ['Animal', /^(animal|mouse|subject|id\b|tumou?r[ _]?id)/i],
  ['Group', /(group|arm|treatment|cohort)/i],
  ['Day', /^(day|days|study[ _]?day|time)\b/i],
  ['Volume', /vol/i],
];
/** A body weight column by its name: any weight that is not the tumor's (ui.ts's isBodyWeight, what the form's
 *  own pick is kept by: the hub must not say "No body weight" of a column Monitor then uses as one). */
export const WEIGHT = /^(?!.*tum).*(weight|\bbw\b)/i;

const role = (name: string) => ROLES.find(([r]) => r === name)![1];

/** A pasted cell as a number: digits with a dot or a comma for the decimals, or with the thousands marked as well
 *  (1,234.5 or 1.234,5). NaN for what cannot be read without guessing: "1,234" alone (1.234 or 1234?), a unit, "<50".
 *  null: an empty cell, or one that says so (NA, n/a, -). */
export function pastedNumber(cell: string): number | null {
  const v = cell.trim();
  if (v === '' || /^(na|n\/a|nan|-|\.)$/i.test(v))
    return null;
  if (/^[+-]?\d+(\.\d+)?$/.test(v))
    return Number(v);
  if (/^[+-]?\d{1,3}(,\d{3})+\.\d+$/.test(v))
    return Number(v.replace(/,/g, ''));
  if (/^[+-]?\d{1,3}(\.\d{3})+,\d+$/.test(v))
    return Number(v.replace(/\./g, '').replace(',', '.'));
  // a decimal comma, unless three digits follow it: that reads as thousands just as well
  return /^[+-]?\d+,\d+$/.test(v) && !/,\d{3}$/.test(v) ? Number(v.replace(',', '.')) : NaN;
}

/** A day written as a number: 0, 7, D7, Day 14, day_21, -1, Day 7 (mm3). Dates are not read. */
export function dayOf(cell: string): number | null {
  const m = /^(?:d|day|days|study day)?[\s_]*(-?\d+(?:\.\d+)?)\s*(?:\([^)]*\))?$/i.exec(cell.trim());
  return m ? Number(m[1]) : null;
}

/** The separator the first lines of a delimited text use most, outside quotes. */
export function separatorOf(text: string): string {
  const head = text.replace(/^\uFEFF/, '').split('\n', 5).join('\n').replace(/"[^"]*"/g, '');
  const count = (c: string) => head.split(c).length - 1;
  return count('\t') > 0 ? '\t' : count(';') > count(',') ? ';' : ',';
}

/** A cell as a number (pastedNumber); where semicolons separate the cells, a comma in a number is its decimals, so
 *  "100,512" is 100.512 there (and refused in a comma file, where it reads as thousands just as well). */
export function numberOf(cell: string, decimalComma: boolean): number | null {
  return decimalComma && /^[+-]?\d+,\d+$/.test(cell.trim()) ? Number(cell.trim().replace(',', '.')) :
    pastedNumber(cell);
}

/** The cells of a delimited text (tabs, commas or semicolons, quotes as a spreadsheet writes them), trimmed. */
export function cellsOf(text: string): string[][] {
  const sep = separatorOf(text);
  const rows: string[][] = [];
  let row: string[] = [];
  let cell = '';
  let quoted = false;
  for (let i = text.startsWith('\uFEFF') ? 1 : 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"')
        cell += text[i++];
      else if (c === '"')
        quoted = false;
      else
        cell += c;
    } else if (c === '"' && cell.trim() === '')
      quoted = true;
    else if (c === sep || c === '\n') {
      row.push(cell.trim());
      cell = '';
      if (c === '\n') {
        rows.push(row);
        row = [];
      }
    } else if (c !== '\r')
      cell += c;
  }
  if (cell.trim() !== '' || row.length > 0)
    rows.push([...row, cell.trim()]);
  return rows;
}

export interface WideRow {animal: string, group: string, day: number, volume: number,
  /** The cell the volume came from, [row, column]: the check that nothing was lost (the study table drops it). */
  from: [number, number]}
export interface WideRead {
  shape: 'by animal' | 'by day',
  rows: WideRow[],
  /** Animals read per group, in the table's order. */
  groups: [string, number][],
  /** Numbers left out on purpose, each with its reason (never a group's name). */
  ignored: string[],
  /** Why nothing was read, with the two shapes that are read; null when read. */
  error: string | null,
}

/** A spreadsheet column's letter (0: A). */
export const letter = (c: number): string => (c >= 26 ? letter(Math.floor(c / 26) - 1) : '') +
  String.fromCharCode(65 + c % 26);
const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? '' : 's'}`;
/** A cell that holds a number, readable or not. */
const isValue = (s: string) => /^[+-]?\d[\d.,]*$/.test(s);
/** A cell that names a role by its word (a name with a number in it is a value, "Mouse 1"). */
const names = (s: string, ...roles: string[]) => !/\d/.test(s) && roles.some((r) => role(r).test(s));
/** A row of the CRO's own statistics, by its label: left out, never an animal. */
const SUMMARY = /^(mean|average|sem|se|sd|std|total|n|median)(\b|=|$)/i;

/** The two shapes read, each as a tiny example, after what was not understood. */
export const SHAPES = 'PMX Studio reads two wide shapes. A column per animal: ",Vehicle,,Cpd-A" / "Day,V1,V2,C1" / ' +
  '"0,37.9,41.2,39.5" (each group\'s name over its first animal). A column per day: "Animal,Group,Day 0,Day 3" / ' +
  '"V1,Vehicle,37.9,51.8" / "V2,,41.2,55.0". Or open it as one row per measurement (animal, group, day, volume).';

/** The plain-words account of a read: animals per group and the days, then what was left out and why. Blinded, the
 *  counts alone, in the table's order. */
export function wideSummary(w: WideRead, blinded = false): string[] {
  if (w.error)
    return [w.error];
  const counts = w.groups.map(([, n]) => String(n));
  const groups = blinded ? `${counts.length > 1 ? `${counts.slice(0, -1).join(', ')} and ` : ''}` +
    `${counts[counts.length - 1] ?? 0} animals` : w.groups.map(([g, n]) => `${g} ${plural(n, 'animal')}`).join(', ');
  const days = w.rows.map((r) => r.day);
  return [`Read ${plural(w.groups.length, 'group')}: ${groups}, days ${Math.min(...days)}–${Math.max(...days)}.`,
    ...w.ignored.map((l) => `Not read: ${l}.`)];
}

/** Every number outside the table's labels (`label`) is read from exactly one cell or left out on purpose
 *  (`ignored`); else why nothing may be read. */
export function unplaced(cells: string[][], rows: WideRow[], label: (r: number, c: number) => boolean,
  ignored: (r: number, c: number) => boolean): string | null {
  const read = new Set(rows.map((x) => `${x.from[0]},${x.from[1]}`));
  const pairs = new Set(rows.map((x) => `${x.animal}\u0000${x.day}`));
  const lost: string[] = [];
  cells.forEach((row, r) => row.forEach((s, c) => {
    if (isValue(s) && !label(r, c) && !ignored(r, c) && !read.has(`${r},${c}`))
      lost.push(`${letter(c)}${r + 1}`);
  }));
  if (pairs.size < rows.length)
    return 'Nothing read: an animal was read twice on one day';
  return lost.length === 0 ? null : `Nothing read: ${plural(lost.length, 'number')} could not be placed ` +
    `(${lost.slice(0, 6).join(', ')}${lost.length > 6 ? ', …' : ''})`;
}

/** A wide table of tumor volumes read into long rows (one per animal and day), in one of two shapes only, or a
 *  refusal that says what was not understood. By animal: a row of group names (each over its group's first animal,
 *  blanks to its right), a row of animal ids headed Day or blank, then a row per day. By day: a header row with the
 *  animal, the group (blank under its name, as merged cells) and the days, then a row per animal. Rows of means or
 *  counts are left out and listed. null when the table does not look wide (the long reading takes it). */
export function readWide(cells: string[][], decimalComma = false): WideRead | null {
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const filled = cells.map((_, r) => r).filter((r) => cells[r].some((s) => s !== ''));
  if (filled.length < 2)
    return null;
  const top = cells[filled[0]];
  // a long table names the day and the animal, the volume or the group
  if (top.some((s) => names(s, 'Day')) && top.some((s) => names(s, 'Animal', 'Volume', 'Group')))
    return null;
  // wide by its first rows: a header of days, or days down column A under a header row (then read strictly)
  const head = filled.slice(0, 4);
  const daysAcross = (r: number) => {
    const days = cells[r].filter((s, c) => c > 0 && s && dayOf(s) !== null);
    return days.length >= 2 && (days.some((s) => !isValue(s)) || cells[r].some((s) => names(s, 'Animal')));
  };
  const dayRow = (r: number) => dayOf(at(r, 0)) !== null && cells[r].some((s, c) => c > 0 && isValue(s));
  const titled = (r: number) => (at(r, 0) === '' || names(at(r, 0), 'Day')) && cells[r].some((s, c) => c > 0 && s);
  if (head.slice(0, 3).some(daysAcross)) {
    return daysAcross(filled[0]) ? readByDay(cells, filled, decimalComma) :
      refuse('by day', 'there are rows above the header row of days');
  }
  return head.some((r, k) => titled(r) && head.slice(k + 1).some(dayRow)) ?
    readByAnimal(cells, filled, decimalComma) : null;
}

const refuse = (shape: WideRead['shape'], what: string): WideRead =>
  ({shape, rows: [], groups: [], ignored: [], error: `Not read: ${what}. ${SHAPES}`});

/** The values of a data row from column `from` on: numbers or blanks; else what is not a number. */
function valuesOf(cells: string[][], r: number, cols: number[], decimalComma: boolean): (number | null)[] | string {
  const out: (number | null)[] = [];
  for (const c of cols) {
    const s = cells[r]?.[c] ?? '';
    const x = numberOf(s, decimalComma);
    if (x !== null && Number.isNaN(x))
      return `cell ${letter(c)}${r + 1} ("${s}") is not a plain number`;
    out.push(x);
  }
  return out;
}

/** The animals' names: their ids, with the group when ids repeat between groups; null when an id repeats inside a
 *  group. */
function animalNames(ids: string[], groups: string[]): string[] | null {
  const inGroup = new Set(ids.map((id, k) => `${groups[k]}\u0000${id}`));
  if (inGroup.size < ids.length)
    return null;
  return new Set(ids).size < ids.length ? ids.map((id, k) => `${groups[k]} ${id}`) : ids;
}

function counted(rows: WideRow[]): [string, number][] {
  const seen = new Map<string, Set<string>>();
  rows.forEach((x) => seen.set(x.group, (seen.get(x.group) ?? new Set()).add(x.animal)));
  return [...seen].map(([g, a]) => [g, a.size]);
}

/** A row of group names, a row of animal ids, then a row per day. */
function readByAnimal(cells: string[][], filled: number[], decimalComma: boolean): WideRead {
  const no = (what: string) => refuse('by animal', what);
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const [g, i, ...rest] = filled;
  if (at(g, 0) !== '' || dayOf(at(i, 0)) !== null)
    return no('the file does not start with a row of group names over a row of animal ids');
  if (at(i, 0) !== '' && !names(at(i, 0), 'Day'))
    return no(`the row of animal ids starts with "${at(i, 0)}", not "Day" or a blank`);
  let last = cells[i].length - 1;
  while (last > 0 && at(i, last) === '')
    last--;
  const cols = Array.from({length: last}, (_, k) => k + 1);
  const missing = cols.find((c) => at(i, c) === '');
  if (cols.length === 0 || missing !== undefined)
    return no(`column ${letter(missing ?? 1)} has no animal id`);
  if (at(g, 1) === '')
    return no('the first animal (column B) has no group name above it');
  const beyond = cells[g].findIndex((s, c) => c > last && s !== '');
  if (beyond >= 0)
    return no(`the group name "${at(g, beyond)}" (column ${letter(beyond)}) is over no animal`);
  const groupOf: string[] = [];
  cols.forEach((c) => groupOf.push(at(g, c) || groupOf[groupOf.length - 1]));
  const measure = groupOf.find((x) => WEIGHT.test(x) || /vol|tum/i.test(x));
  if (measure)
    return no(`"${measure}" over the animals is a measurement, not a group`);
  const ids = cols.map((c) => at(i, c));
  const stat = cols.find((c) => SUMMARY.test(at(i, c)));
  if (stat !== undefined)
    return no(`column ${letter(stat)} is headed "${at(i, stat)}", a summary, not an animal`);
  const animals = animalNames(ids, groupOf);
  if (!animals)
    return no('an animal id is there twice in one group');

  const rows: WideRow[] = [];
  const summary: number[] = [];
  let previous = -Infinity;
  for (const r of rest) {
    if (SUMMARY.test(at(r, 0))) {
      summary.push(r);
      continue;
    }
    if (summary.length > 0)
      return no(`row ${r + 1} comes after the summary rows`);
    const day = dayOf(at(r, 0));
    if (day === null) {
      return no(at(r, 0) === '' ? `row ${r + 1} has no day in column A` :
        `row ${r + 1} starts with "${at(r, 0)}", not a day (0, 7, D7 or Day 7; dates are not read)`);
    }
    if (day <= previous)
      return no(`day ${at(r, 0)} (row ${r + 1}) does not come after the day above it`);
    previous = day;
    const beyondCell = cells[r].findIndex((s, c) => c > last && s !== '');
    if (beyondCell >= 0)
      return no(`cell ${letter(beyondCell)}${r + 1} has no animal over it`);
    const values = valuesOf(cells, r, cols, decimalComma);
    if (typeof values === 'string')
      return no(values);
    values.forEach((v, k) => v !== null &&
      rows.push({animal: animals[k], group: groupOf[k], day, volume: v, from: [r, cols[k]]}));
  }
  if (rows.length === 0)
    return no('no measurements under the animals');
  const left = unplaced(cells, rows, (r, c) => r === g || r === i || c === 0, (r) => summary.includes(r));
  return {shape: 'by animal', rows, groups: counted(rows), error: left,
    ignored: summary.map((r) => `row ${r + 1} ("${at(r, 0)}"), a summary row`)};
}

/** A header row (animal, group, days), then a row per animal. */
function readByDay(cells: string[][], filled: number[], decimalComma: boolean): WideRead {
  const no = (what: string) => refuse('by day', what);
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const [h, ...rest] = filled;
  const head = cells[h];
  const first = head.findIndex((s) => s !== '' && dayOf(s) !== null);
  const labels = Array.from({length: first}, (_, c) => c);
  const animalCols = labels.filter((c) => names(head[c], 'Animal'));
  const groupCols = labels.filter((c) => names(head[c], 'Group') && !animalCols.includes(c));
  if (animalCols.length !== 1)
    return no('no column is headed as the animals (Animal, Mouse, ID)');
  if (groupCols.length !== 1) {
    return no(groupCols.length === 0 ? 'no column is headed as the groups (Group, Arm, Treatment)' :
      `more than one column is headed as the groups (${groupCols.map((c) => `"${head[c]}"`).join(', ')})`);
  }
  const [animalCol] = animalCols;
  const [groupCol] = groupCols;
  // the days, all of them going up; blank-headed columns after them only
  let lastDay = first;
  for (let c = first; c < head.length; c++) {
    if (head[c] === '')
      continue;
    if (dayOf(head[c]) === null)
      return no(`column ${letter(c)} is headed "${head[c]}", not a day (0, 7, D7 or Day 7)`);
    if (c > first && dayOf(head[c])! <= dayOf(head[lastDay])!)
      return no(`the day "${head[c]}" (column ${letter(c)}) does not come after the one before it`);
    if (c > lastDay + 1)
      return no(`column ${letter(lastDay + 1)} between the days has no heading`);
    lastDay = c;
  }
  const dayCols = Array.from({length: lastDay - first + 1}, (_, k) => first + k);
  const width = Math.max(...filled.map((r) => cells[r].length));
  const after = Array.from({length: Math.max(0, width - lastDay - 1)}, (_, k) => lastDay + 1 + k)
    .filter((c) => rest.some((r) => isValue(at(r, c))));

  const rows: WideRow[] = [];
  const summary: number[] = [];
  const animals: {id: string, group: string, r: number}[] = [];
  let group = '';
  for (const r of rest) {
    const aside = labels.find((c) => c !== animalCol && c !== groupCol && isValue(at(r, c)));
    if (SUMMARY.test(at(r, animalCol)) || SUMMARY.test(at(r, groupCol)) || at(r, animalCol) === '' &&
      at(r, groupCol) === '' && SUMMARY.test(cells[r].find((s) => s) ?? '')) {
      summary.push(r);
      continue;
    }
    if (aside !== undefined)
      return no(`column ${letter(aside)} ("${head[aside]}") holds numbers before the days`);
    if (at(r, groupCol) !== '')
      group = at(r, groupCol);
    if (at(r, animalCol) === '')
      return no(`row ${r + 1} has no animal id`);
    if (group === '')
      return no(`row ${r + 1} ("${at(r, animalCol)}") has no group`);
    const values = valuesOf(cells, r, dayCols, decimalComma);
    if (typeof values === 'string')
      return no(values);
    if (values.every((v) => v === null))
      return no(`row ${r + 1} ("${at(r, animalCol)}") has no measurements`);
    animals.push({id: at(r, animalCol), group, r});
  }
  const named = animalNames(animals.map((a) => a.id), animals.map((a) => a.group));
  if (!named)
    return no('an animal id is there twice in one group');
  animals.forEach((a, k) => {
    const values = valuesOf(cells, a.r, dayCols, decimalComma) as (number | null)[];
    values.forEach((v, j) => v !== null && rows.push({animal: named[k], group: a.group, day: dayOf(head[dayCols[j]])!,
      volume: v, from: [a.r, dayCols[j]]}));
  });
  if (rows.length === 0)
    return no('no measurements under the days');
  const left = unplaced(cells, rows, (r, c) => r === h || c < first, (r, c) => summary.includes(r) ||
    after.includes(c));
  return {shape: 'by day', rows, groups: counted(rows), error: left, ignored: [
    ...summary.map((r) => `row ${r + 1} ("${cells[r].find((s) => s)}"), a summary row`),
    ...after.map((c) => `column ${letter(c)}, with no heading, after the days`)]};
}
