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

/** A day written as a number: 0, 7, D7, Day 14, day_21, -1. Dates are not read. */
export function dayOf(cell: string): number | null {
  const m = /^(?:d|day|days|study day)?[\s_]*(-?\d+(?:\.\d+)?)$/i.exec(cell.trim());
  return m ? Number(m[1]) : null;
}

/** The cells of a delimited text (tabs, commas or semicolons, quotes as a spreadsheet writes them), trimmed. The
 *  separator is the one the first lines use most, outside quotes. */
export function cellsOf(text: string): string[][] {
  const head = text.replace(/^\uFEFF/, '').split('\n', 5).join('\n').replace(/"[^"]*"/g, '');
  const count = (c: string) => head.split(c).length - 1;
  const sep = count('\t') > 0 ? '\t' : count(';') > count(',') ? ';' : ',';
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

export interface WideRow {animal: string, group: string, day: number, volume: number | null, weight: number | null}
/** What only the user can say: `choices` by their index, given back in `Answers[id]`; `blind`, the same asked without
 *  naming a group, for a blinded study. */
export interface Question {id: string, text: string, choices: string[], blind?: {text: string, choices: string[]}}
export type Answers = {[id: string]: number};
export interface WideRead {
  shape: 'by day' | 'by animal',
  rows: WideRow[],
  /** Animals read per group, in the table's order. */
  groups: [string, number][],
  /** Every group name the table gives, read or not. */
  named: string[],
  /** Measurements not read, by where they are and why (never by a group's name). */
  left: string[],
  /** How the table was read, where it needed telling. */
  notes: string[],
  question: Question | null,
  error: string | null,
}

/** A spreadsheet column's letter (0: A). */
export const letter = (c: number): string => (c >= 26 ? letter(Math.floor(c / 26) - 1) : '') +
  String.fromCharCode(65 + c % 26);
const plural = (n: number, word: string) => `${n} ${word}${n === 1 ? '' : 's'}`;
const list = (xs: string[], max = 6) => xs.slice(0, max).join(', ') + (xs.length > max ? `, … (${xs.length})` : '');
/** A cell that holds a value, readable or not. */
const isValue = (s: string) => /^[+-]?\d[\d.,]*$/.test(s);
/** A cell that names a role by its word (a name with a number in it is a value, "Mouse 1"). */
const names = (s: string, ...roles: string[]) => !/\d/.test(s) && roles.some((r) => role(r).test(s));
const NO_NAME = 'No group name';
const rowsOf = (rs: number[]) => rs.length === 1 ? `row ${rs[0] + 1}` : `rows ${rs[0] + 1}–${rs[rs.length - 1] + 1}`;

/** The plain-words account of a read: animals per group first. Blinded, no group is named: the counts come in the
 *  table's order. */
export function wideSummary(w: WideRead, blinded = false): string[] {
  if (w.error)
    return [w.error];
  const counts = w.groups.map(([, n]) => String(n));
  const groups = blinded ? `${counts.length > 1 ? `${counts.slice(0, -1).join(', ')} and ` : ''}` +
    `${counts[counts.length - 1] ?? 0} animals` : w.groups.map(([g, n]) => `${g} ${plural(n, 'animal')}`).join(', ');
  const notes = blinded ? w.notes.filter((n) => !w.named.some((g) => n.includes(`"${g}"`))) : w.notes;
  return [`Read ${plural(w.groups.length, 'group')}: ${groups}.`, ...w.left.map((l) => `Not read: ${l}.`), ...notes];
}

/** A wide table of tumor volumes, as most CROs send it, read into long rows (one per animal and day): by day (a row
 *  per animal, a column per day, the group in a column or in a row over its animals) or by animal (a row per day, a
 *  column per animal, each group's name over its first animal only, as merged cells come out of a spreadsheet).
 *  Every measurement ends up in a group or is named in what is not read; where the groups cannot be told without
 *  guessing, a question comes back instead of rows; a refusal says why. null when the cells are not a wide table. */
export function readWide(cells: string[][], answers: Answers = {}): WideRead | null {
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const width = Math.max(0, ...cells.map((r) => r.length));
  const valueCols = (r: number, from = 0) => {
    const out: number[] = [];
    for (let c = from; c < width; c++) {
      if (isValue(at(r, c)))
        out.push(c);
    }
    return out;
  };
  const firstData = cells.findIndex((_, r) => valueCols(r).length > 0);
  // a long table names the day and the animal, the volume or the group; a summary gives means
  const head = cells.slice(0, Math.max(1, firstData)).flat();
  const top = cells[0] ?? [];
  if (firstData < 0 || top.some((s) => names(s, 'Day')) && top.some((s) => names(s, 'Animal', 'Volume', 'Group')) ||
    head.some((s) => /^(mean|average|sem|se|sd|std|s\.e\.m\.?)\b/i.test(s)))
    return null;

  // by day: a header row whose last cells are days, after labels that are not numbers, with nothing measured above
  let byDay: {r: number, k: number} | null = null;
  for (let r = 0; r <= firstData && !byDay; r++) {
    const filled = cells[r].map((s, c) => [s, c] as [string, number]).filter(([s]) => s !== '');
    let k = filled.length;
    while (k > 0 && dayOf(filled[k - 1][0]) !== null)
      k--;
    const labels = filled.slice(0, k).map(([s]) => s);
    // days written as plain numbers look like a row of measurements, unless labelled as animals or groups
    const labelled = r < firstData || labels.some((s) => names(s, 'Animal', 'Group'));
    if (filled.length - k >= 2 && filled[k][1] >= 1 && labels.every((s) => !isValue(s)) &&
      (labelled || labels.length === 0))
      byDay = {r, k: filled[k][1]};
  }
  // by animal: the days down the leftmost column, or down a column headed Day, under header rows that put no id or
  // group name over them; a day measured twice is named, a day going back ends the run
  let byAnimal: ByAnimal | null = null;
  const lead = Math.min(...cells.map((row) => row.findIndex((s) => s !== '')).filter((c) => c >= 0));
  const headed = Array.from({length: width}, (_, c) => c)
    .filter((c) => cells.some((row) => names(row[c] ?? '', 'Day')));
  for (const c of [lead, ...headed]) {
    const start = cells.findIndex((_, r) => dayOf(at(r, c)) !== null && valueCols(r, c + 1).length > 0);
    if (byAnimal || start < 1 || !cells.slice(0, start).every((row) => !row[c] || names(row[c], 'Day')) ||
      cells[start - 1].some((s) => names(s, 'Animal', 'Volume')))
      continue;
    const rows: number[] = [];
    const after: number[] = [];
    const again: number[] = [];
    for (let r = start; r < cells.length; r++) {
      if (!cells[r].some((s, k) => k > c && s !== ''))
        continue;
      const d = dayOf(at(r, c));
      const last = rows.length ? dayOf(at(rows[rows.length - 1], c))! : -Infinity;
      if (after.length === 0 && d !== null && d > last)
        rows.push(r);
      else if (after.length === 0 && d === last)
        again.push(r);
      else
        after.push(r);
    }
    // most rows in one run of days going up: a long table restarts its days with each animal
    if (rows.length >= 2 && rows.length > after.length)
      byAnimal = {c, rows, after, again};
  }
  if (!byDay && !byAnimal)
    return null;
  if (byDay && byAnimal) {
    const labels = cells[byDay.r].slice(0, byDay.k);
    const dayHead = names(at(byAnimal.rows[0] - 1, byAnimal.c), 'Day') || labels.some((s) => names(s, 'Day'));
    const animalHead = labels.some((s) => names(s, 'Animal', 'Group'));
    if (dayHead === animalHead) {
      const shape = answers.shape;
      if (shape === undefined) {
        return ask('by day', {id: 'shape', text: 'Are the days down the first column (one row per day, a column ' +
          'per animal) or across the top (one row per animal, a column per day)?',
        choices: ['Down the first column: one row per day', 'Across the top: one row per animal']});
      }
      return shape === 0 ? readByAnimal(cells, byAnimal, answers) : readByDay(cells, byDay, answers);
    }
    return dayHead ? readByAnimal(cells, byAnimal, answers) : readByDay(cells, byDay, answers);
  }
  return byDay ? readByDay(cells, byDay, answers) : readByAnimal(cells, byAnimal!, answers);
}

interface ByAnimal {c: number, rows: number[], after: number[], again: number[]}

function ask(shape: WideRead['shape'], question: Question): WideRead {
  return {shape, rows: [], groups: [], named: [], left: [], notes: [], question, error: null};
}
function refuse(shape: WideRead['shape'], error: string): WideRead {
  return {shape, rows: [], groups: [], named: [], left: [], notes: [], question: null, error};
}

/** The question for measurements that come before the first group name, or the answer to it: the group they go to
 *  (null: left out). */
function beforeFirst(where: string, first: string, answers: Answers): {question?: Question, group: string | null} {
  const a = answers.before;
  if (a === undefined) {
    const text = `${where} hold measurements but come before the first group name. Which group are they in?`;
    const rest = [`A group of their own ("${NO_NAME}")`, 'None: leave them out'];
    return {group: null, question: {id: 'before', text, choices: [`${first}, the group after them`, ...rest],
      blind: {text, choices: ['The group after them', ...rest]}}};
  }
  return {group: a === 0 ? first : a === 1 ? NO_NAME : null};
}

/** The rows of a read: the values, the numbers that cannot be read (a refusal), the words that are not values, and
 *  the animals that brought none (named where they are). */
class Collect {
  rows = new Map<string, WideRow>();
  read = new Set<string>();
  unread: string[] = [];
  words: string[] = [];
  counts = new Map<string, number>();

  add(animal: string, group: string, day: number, r: number, c: number, cell: string, what: 'volume' | 'weight'):
    void {
    const x = pastedNumber(cell);
    if (cell !== '' && x !== null && Number.isNaN(x))
      (/\d/.test(cell) ? this.unread : this.words).push(`"${cell}" at ${letter(c)}${r + 1}`);
    if (x === null || Number.isNaN(x))
      return;
    const key = `${animal}\u0000${day}`;
    if (!this.read.has(animal)) {
      this.read.add(animal);
      this.counts.set(group, (this.counts.get(group) ?? 0) + 1);
    }
    const row = this.rows.get(key) ?? {animal, group, day, volume: null, weight: null};
    row[what] = x;
    this.rows.set(key, row);
  }

  /** `animals`: each one's name, the id the file gives it ('' for none), and where it is. */
  done(shape: WideRead['shape'], named: string[], left: string[], notes: string[],
    animals: {name: string, id: string, where: string}[]): WideRead {
    if (this.unread.length > 0) {
      return refuse(shape, `Nothing read: ${plural(this.unread.length, 'value')} could not be read as numbers ` +
        `(${list(this.unread, 4)}). Write them as plain numbers, without units or thousands separators`);
    }
    if (this.words.length > 0)
      left.push(`cells that are not numbers: ${list(this.words, 4)}`);
    for (const a of animals.filter((x) => !this.read.has(x.name)))
      left.push(`${a.where}${a.id ? ` (${a.id})` : ''}, no numbers`);
    return {shape, rows: [...this.rows.values()], groups: [...this.counts], named, left, notes, question: null,
      error: null};
  }
}

/** Names the animals: their ids, with the group when ids repeat between groups; a refusal when one repeats inside a
 *  group (two blocks side by side, most often). */
function nameAnimals(animals: {id: string, group: string, where: string}[], notes: string[]): string[] | string {
  const seen = new Map<string, string>();
  for (const a of animals) {
    const key = `${a.group}\u0000${a.id}`;
    if (seen.has(key)) {
      return `Nothing read: animal "${a.id}" is there twice in one group (${seen.get(key)} and ${a.where}). If ` +
        'the table holds body weight or another measurement beside the volumes, open the volumes on their own';
    }
    seen.set(key, a.where);
  }
  const groupsOf = new Map<string, Set<string>>();
  animals.forEach((a) => groupsOf.set(a.id, (groupsOf.get(a.id) ?? new Set()).add(a.group)));
  if ([...groupsOf.values()].some((g) => g.size > 1)) {
    notes.push('Animal ids repeat between groups: each animal is named with its group.');
    return animals.map((a) => `${a.group} ${a.id}`);
  }
  return animals.map((a) => a.id);
}

/** One row per day, a column per animal. */
function readByAnimal(cells: string[][], p: ByAnimal, answers: Answers): WideRead {
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const width = Math.max(0, ...cells.map((r) => r.length));
  const first = p.rows[0];
  const left: string[] = [];
  const notes: string[] = [];
  // an animal's column: anything under the header in the run of days, a number or not
  const measured: number[] = [];
  for (let c = 0; c < width; c++) {
    if (c < p.c && p.rows.some((r) => isValue(at(r, c))))
      left.push(`column ${letter(c)}${at(first - 1, c) ? ` ("${at(first - 1, c)}")` : ''}, left of the days`);
    else if (c > p.c && p.rows.some((r) => at(r, c) !== ''))
      measured.push(c);
  }
  for (const r of p.again)
    left.push(`row ${r + 1}, day ${at(r, p.c)} a second time`);
  if (p.after.length > 0) {
    left.push(`${rowsOf(p.after)}, where the days stop going up (row ` +
      `${p.after[0] + 1} starts "${at(p.after[0], p.c)}")`);
  }
  // the row over the measurements: animal ids when a row above it names something, and it misses few ids
  const idRow = first - 1;
  const above = (r: number) => (cells[r] ?? []).some((s, c) => c > p.c && s !== '');
  let groupRow = -1;
  for (let r = idRow - 1; r >= 0 && groupRow < 0; r--) {
    if (above(r))
      groupRow = r;
  }
  const blanks = measured.filter((c) => at(idRow, c) === '');
  let ids = groupRow >= 0 && blanks.length <= measured.length / 2;
  if (!ids) {
    groupRow = idRow;
    if (blanks.length === 0 && new Set(measured.map((c) => at(idRow, c))).size === measured.length) {
      // every column named, each name once: animals of one group, or groups of one animal each
      const a = answers.header;
      if (a === undefined) {
        const choices = ['The animals, all in one group', 'The groups, one animal each'];
        return ask('by animal', {id: 'header', text: `Is the row over the measurements ("${at(idRow, measured[0])}"` +
          `, "${at(idRow, measured[1])}", …) the animals or their groups?`, choices,
        blind: {text: 'Is the row over the measurements the animals or their groups?', choices}});
      }
      if (a === 0) {
        ids = true;
        groupRow = -1;
      }
    }
  }
  // each column's group: the last name at or left of it (merged cells come through as the name, then blanks)
  const groupOf = new Map<number, string>();
  let current = '';
  for (let c = p.c + 1; c < width; c++) {
    if (groupRow >= 0 && at(groupRow, c) !== '')
      current = at(groupRow, c);
    groupOf.set(c, groupRow < 0 ? 'One group' : current);
  }
  const named = [...new Set(groupOf.values())].filter((g) => g !== '');
  const measures = named.filter((g) => WEIGHT.test(g) || /vol|tum/i.test(g));
  if (measures.length > 0) {
    return refuse('by animal', `Nothing read: the names over the animals (${list(measures.map((g) => `"${g}"`))}) ` +
      'are measurements, not groups. Open the tumor volumes on their own, each group\'s name over its first animal');
  }
  const before = measured.filter((c) => groupOf.get(c) === '');
  if (before.length > 0) {
    const cols = before.map(letter);
    const b = beforeFirst(`Column${cols.length === 1 ? '' : 's'} ${list(cols)}`, named[0] ?? NO_NAME, answers);
    if (b.question)
      return ask('by animal', b.question);
    before.forEach((c) => b.group === null ? groupOf.delete(c) : groupOf.set(c, b.group));
    if (b.group === null)
      left.push(`column${cols.length === 1 ? '' : 's'} ${list(cols)}, before the first group name (as you chose)`);
  }
  const kept = measured.filter((c) => groupOf.has(c));
  const numbered = new Map<string, number>();
  const raw = kept.map((c) => {
    const group = groupOf.get(c)!;
    numbered.set(group, (numbered.get(group) ?? 0) + 1);
    const id = ids ? at(idRow, c) : '';
    return {id: id || (ids ? `column ${letter(c)}` : `${group} ${numbered.get(group)}`), group, file: id,
      where: `column ${letter(c)}`};
  });
  if (!ids)
    notes.push('The animals have no ids in the file: they are numbered within each group.');
  else if (blanks.length > 0) {
    notes.push(`No animal id over column${blanks.length === 1 ? '' : 's'} ${list(blanks.map(letter))}: named by ` +
      'the column.');
  }
  const animals = nameAnimals(raw, notes);
  if (typeof animals === 'string')
    return refuse('by animal', animals);
  const empty = named.filter((g) => !kept.some((c) => groupOf.get(c) === g) && g !== NO_NAME);
  if (empty.length > 0)
    notes.push(`No measurements under ${list(empty.map((g) => `"${g}"`))}.`);
  const out = new Collect();
  kept.forEach((c, k) => p.rows.forEach((r) => out.add(animals[k], groupOf.get(c)!, dayOf(at(r, p.c))!, r, c,
    at(r, c), 'volume')));
  notes.push('Body weight is not read from a table with a column per animal.');
  return out.done('by animal', named, left, notes, raw.map((a, k) => ({name: animals[k], id: a.file,
    where: a.where})));
}

/** One row per animal, a column per day; the group in a column of its own, or in a row over its animals. */
function readByDay(cells: string[][], p: {r: number, k: number}, answers: Answers): WideRead {
  const at = (r: number, c: number) => cells[r]?.[c] ?? '';
  const width = Math.max(0, ...cells.map((r) => r.length));
  const left: string[] = [];
  const notes: string[] = [];
  // the days, in blocks: a block ends where the days start again
  const blocks: {cols: number[], title: string}[] = [];
  for (let c = p.k; c < width; c++) {
    const d = dayOf(at(p.r, c));
    if (d === null)
      continue;
    const last = blocks[blocks.length - 1];
    if (!last || d <= dayOf(at(p.r, last.cols[last.cols.length - 1]))!)
      blocks.push({cols: [c], title: ''});
    else
      last.cols.push(c);
  }
  for (const b of blocks) {
    for (let r = p.r - 1; r >= 0 && !b.title; r--)
      b.title = b.cols.map((c) => at(r, c)).find((s) => s !== '') ?? '';
  }
  // the rows of animals: down to a second header row, if any (another table below)
  let end = cells.length;
  for (let r = p.r + 1; r < cells.length && end === cells.length; r++) {
    const again = cells[r].slice(p.k).filter((s) => s !== '' && !isValue(s) && dayOf(s) !== null).length >= 2 ||
      cells[r].slice(0, p.k).some((s, c) => s !== '' && s === at(p.r, c));
    if (again)
      end = r;
  }
  const valued = (r: number) => cells[r].some((s, c) => c >= p.k && isValue(s));
  const below = Array.from({length: cells.length - end}, (_, k) => end + k).filter(valued);
  if (below.length > 0)
    left.push(`${rowsOf(below)}, a second table below the first`);
  const rows = Array.from({length: end - p.r - 1}, (_, k) => p.r + 1 + k).filter((r) => cells[r].some((s) => s));
  const dayCols = new Set(blocks.flatMap((b) => b.cols));
  for (let c = p.k; c < width; c++) {
    if (!dayCols.has(c) && rows.some((r) => isValue(at(r, c))))
      left.push(`column ${letter(c)}, which has no day above it`);
  }

  // which block is the volume, which the body weight
  const span = (b: {cols: number[], title: string}) => `columns ${letter(b.cols[0])}–` +
    `${letter(b.cols[b.cols.length - 1])}${b.title ? ` ("${b.title}")` : ''}`;
  const isWeight = (b: {title: string}) => WEIGHT.test(b.title);
  const vols = blocks.filter((b) => !isWeight(b) && /vol|tum/i.test(b.title));
  const weights = blocks.filter(isWeight);
  let volume = vols.length === 1 ? vols[0] : blocks.length === 1 ? blocks[0] :
    vols.length === 0 && weights.length === 1 && blocks.length === 2 ? blocks.find((b) => !isWeight(b))! : null;
  if (volume && isWeight(volume))
    return refuse('by day', `Nothing read: the days hold body weight ("${volume.title}"), not tumor volumes`);
  if (!volume) {
    if (answers.block === undefined) {
      return ask('by day', {id: 'block', text: 'The days come in more than one block. Which one holds the tumor ' +
        'volumes?', choices: blocks.map(span)});
    }
    volume = blocks[answers.block];
  }
  const weight = weights.length === 1 && weights[0] !== volume ? weights[0] : null;
  for (const b of blocks) {
    if (b !== volume && b !== weight)
      left.push(`${span(b)}, not clearly tumor volume or body weight`);
  }

  // the labels left of the days: the animal, the group (a merged cell comes through as the name, then blanks)
  const labels = Array.from({length: p.k}, (_, c) => c);
  const measuredRows = rows.filter(valued);
  const values = (c: number) => measuredRows.map((r) => at(r, c));
  const distinct = (c: number) => values(c).every((s) => s !== '') && new Set(values(c)).size === measuredRows.length;
  let groupCol = labels.find((c) => names(at(p.r, c), 'Group'));
  const animalCol = labels.find((c) => c !== groupCol && names(at(p.r, c), 'Animal')) ??
    labels.find((c) => c !== groupCol && distinct(c));
  if (groupCol === undefined) {
    // no column is headed as the groups: a merged one without a heading is, any other is asked about
    const could = labels.filter((c) => c !== animalCol && values(c).some((s) => s !== '') && !distinct(c));
    if (could.length === 1 && at(p.r, could[0]) === '')
      groupCol = could[0];
    else if (could.length > 0) {
      if (answers.group === undefined) {
        return ask('by day', {id: 'group', text: 'No column is headed as the groups. Which one holds them?',
          choices: [...could.map((c) => `Column ${letter(c)}${at(p.r, c) ? ` ("${at(p.r, c)}")` : ''}`),
            'None: each group is named in a row over its animals']});
      }
      groupCol = could[answers.group];
    }
  }
  if (animalCol === undefined)
    notes.push('No column of animal ids: the animals are numbered within each group.');

  // a row of its own names a group, alone left of the days and not a number (that is an animal without measurements)
  const title = (r: number) => groupCol === undefined && cells[r].every((s, c) => !s || c < p.k) &&
    !isValue(cells[r].find((s) => s)!);
  const raw: {id: string, group: string, where: string, file: string, r: number}[] = [];
  const before: number[] = [];
  let group = '';
  for (const r of rows) {
    if (title(r)) {
      group = cells[r].find((s) => s)!;
      continue;
    }
    if (groupCol !== undefined && at(r, groupCol) !== '')
      group = at(r, groupCol);
    if (!group)
      before.push(r);
    const id = animalCol !== undefined ? at(r, animalCol) : '';
    raw.push({id, file: id, group, where: `row ${r + 1}`, r});
  }
  if (raw.every((a) => !a.group))
    return refuse('by day', 'Nothing read: no group names found, in a column or in a row over the animals');
  if (before.length > 0) {
    const what = `Row${before.length === 1 ? '' : 's'} ${list(before.map((r) => String(r + 1)))}`;
    const b = beforeFirst(what, raw.find((a) => a.group)!.group, answers);
    if (b.question)
      return ask('by day', b.question);
    for (const a of raw.filter((x) => !x.group))
      a.group = b.group ?? '';
    if (b.group === null)
      left.push(`${what.toLowerCase()}, before the first group name (as you chose)`);
  }
  const kept = raw.filter((a) => a.group);
  const numbered = new Map<string, number>();
  for (const a of kept) {
    numbered.set(a.group, (numbered.get(a.group) ?? 0) + 1);
    if (!a.id)
      a.id = animalCol !== undefined ? `row ${a.r + 1}` : `${a.group} ${numbered.get(a.group)}`;
  }
  const nameless = kept.filter((a) => animalCol !== undefined && !at(a.r, animalCol)).map((a) => String(a.r + 1));
  if (nameless.length > 0)
    notes.push(`No animal id in row${nameless.length === 1 ? '' : 's'} ${list(nameless)}: named by the row.`);
  const animals = nameAnimals(kept, notes);
  if (typeof animals === 'string')
    return refuse('by day', animals);
  const out = new Collect();
  kept.forEach((a, k) => {
    for (const c of volume!.cols)
      out.add(animals[k], a.group, dayOf(at(p.r, c))!, a.r, c, at(a.r, c), 'volume');
    for (const c of weight?.cols ?? [])
      out.add(animals[k], a.group, dayOf(at(p.r, c))!, a.r, c, at(a.r, c), 'weight');
  });
  if (!weight && blocks.length === 1)
    notes.push('No body weight in the file.');
  const named = [...new Set(raw.map((a) => a.group).filter((g) => g))];
  return out.done('by day', named, left, notes, kept.map((a, k) => ({name: animals[k], id: a.file,
    where: a.where})));
}
