import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {EntryAnimal, enterDayDialog, entrySummary, loadDraft, pastedRows, readDay, removeDay, rowChecks, saveDraft,
  stagedVolume, summaryText,
  valuesOfDay, writeDay} from '../day-entry';
import {modalOf} from '../ui';
import {checkReading, monitorFrames, outsideBand, removedAnimals, twoReadings} from '../monitor';
import {analysisDeclarationOf} from '../analyze';
import {NOT_A_NUMBER, checkTyped} from '../reading-check';
import {pastedNumber} from '../wide-table';
import {animalMarks, correctValues, lastReadingRow, mark, readingColumns, tagStudy, withoutExcluded}
  from '../study-table';
import {logOf, roleColumns} from '../study-file';
import {emptyMeta, setMetaOn, metaOf} from '../study-meta';
import {exampleStudy} from './example-study';

/** Two animals measured on day 0, a body weight column. */
function study(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['a1', 'a2']),
    DG.Column.fromStrings('group', ['Vehicle', 'Drug']), DG.Column.fromList(DG.TYPE.FLOAT, 'day', [0, 0]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'volume', [100, 120]),
    DG.Column.fromList(DG.TYPE.FLOAT, 'body_weight', [20, 21])]);
  setMetaOn(df, emptyMeta());
  return df;
}
const columnsOf = (df: DG.DataFrame) => {
  const c = roleColumns(df);
  return {animal: c.animal, group: c.group, day: c.day, readings: readingColumns(df, c.volume, df.col('body_weight'))};
};
const animal = (a: string, group: string, last: number): EntryAnimal => ({animal: a, group, lastDay: 0,
  lastVolume: last, lastWeight: 20, lower: 50, upper: 300, expected: '80–250'});

category('Monitor: entering a day', () => {
  test('one rule for a value just measured: ok, measure again, likely typo, impossible', async () => {
    expect(checkReading(150, {before: 100, lower: 50, upper: 300})?.kind, 'ok');
    expect(checkReading(400, {before: 100, lower: 50, upper: 300})?.kind, 'again');
    expect(checkReading(1000, {before: 100, lower: 50, upper: 300})?.kind, 'typo');
    expect(checkReading(-5, {before: 100})?.kind, 'impossible');
    // a tumor gone is a value, not an impossible one; a body weight of 0 is impossible
    expect(checkReading(0, {})?.kind, 'ok');
    expect(checkReading(0, {weight: true})?.kind, 'impossible');
    expect(checkReading(null, {}), null);
    // Monitor's "measure again" list reads the same rule
    expect(outsideBand([{animal: 'a', day: 3, volume: 400}], 3, new Map([['a', {check_lower: 50,
      check_upper: 300}]])).length, 1);
    expect(rowChecks(animal('a1', 'V', 100), {volume: '12,x'}, false).volume?.kind, 'impossible');
    // one reading for an entry and a correction (QA 0.1.3 round 3): a decimal comma reads, "1e9" is no number, a
    // volume bigger than a mouse is impossible, not a typo to check
    expect(checkTyped(pastedNumber('289,3'), {before: 250})?.kind, 'ok');
    expect(checkTyped(pastedNumber('1e9'), {})?.why, NOT_A_NUMBER.why);
    const huge = checkTyped(pastedNumber('1,000,000,000.0'), {before: 250});
    expect([huge?.kind, huge?.why].join(' | '), 'impossible | 1,000,000,000 mm³ is not a possible tumor volume');
    expect(checkTyped(5000, {weight: true})?.kind, 'impossible');
    // (QA 0.1.3 round 4) 50,000 mm³ is no volume at all; 0 after a value asks, never refused; 2,000 g still a typo
    expect(checkTyped(50000, {before: 250})?.kind, 'impossible');
    expect([checkTyped(0, {before: 250})?.kind, checkTyped(0, {before: 250})?.why.startsWith('0 mm³: tumor gone?'),
      checkTyped(0, {before: 0})?.kind].join(), 'again,true,ok');
    expect(checkTyped(2000, {weight: true, before: 20})?.kind, 'typo');
  });

  test('the footer counts what is entered and what to check; calipers give the volume', async () => {
    const list = [animal('a1', 'V', 100), animal('a2', 'D', 120), animal('a3', 'D', 90)];
    const s = entrySummary(list, {a1: {volume: '150'}, a2: {volume: '1000'}, a3: {missed: 'Dead'}}, false);
    expect(summaryText(s, 3), '3 of 3 entered · 1 to check');
    expect(stagedVolume({length: '10', width: '6'}, true), 180);
    expect(stagedVolume({length: '10'}, true), null);
    expect(pastedRows('120\r\n130,5\n\n').map((r) => r[0]).join('|'), '120|130,5');
    // two columns pasted (volume, then body weight): both kept, a row each
    expect(pastedRows('120\t20,1\n130\t19.8\n').map((r) => r.join('/')).join('|'), '120/20,1|130/19.8');
    // an animal off the study is listed, not counted
    const off = {...animal('a4', 'D', 90), gone: 'Off the study: dead on day 2'};
    expect(summaryText(entrySummary([...list, off], {a1: {volume: '150'}}, false),
      [...list, off].filter((a) => !a.gone).length), '1 of 3 entered');
  });

  test('a day typed: whole days from 0, a day far past the study asked first', async () => {
    const o = {last: 9, end: 21};
    expect(readDay('7', o).day, 7);
    expect(readDay(' ', o).error, 'Type the study day (0 = randomization).');
    expect(readDay('2.5', o).error?.startsWith('Whole days only'), true);
    expect(readDay('abc', o).error?.includes('not a day'), true);
    expect(readDay('-3', o).error?.includes('0 or later'), true);
    expect(readDay('9999', o).ask?.startsWith('Day 9999 is far past the study'), true);
    expect(readDay('30', o).ask === undefined, true);
    // a day still to come asks first too (QA 0.1.3 round 4 chaos N9)
    const t = {...o, today: 12};
    expect([readDay('12', t).ask, readDay('13', t).ask?.startsWith('Day 13 is tomorrow (today is day 12)')].join(),
      ',true');
  });

  test('a saved day can be taken back: a value erased, or the whole day, the old values in the log', async () => {
    const df = study();
    const c = columnsOf(df);
    const list = [animal('a1', 'Vehicle', 100), animal('a2', 'Drug', 120)];
    writeDay(df, c, 3, list, {a1: {volume: '150', weight: '20'}, a2: {volume: '170'}});
    expect(df.rowCount, 4);
    // a2's day 3 erased: its row goes, the log says what it held
    writeDay(df, c, 3, list, {a1: {volume: '150', weight: '20'}});
    expect(df.rowCount, 3);
    expect(logOf(df).pop()!.what.includes('1 value removed (a2 170 mm³)'), true);
    // a1's weight alone erased: the row stays, its weight empty
    writeDay(df, c, 3, list, {a1: {volume: '150'}});
    expect(df.col('body_weight')!.isNone(2), true);
    // the whole day
    expect(removeDay(df, c, 3), 'Day 3 removed (1 value) · values: a1 150 mm³');
    expect(df.rowCount, 2);
    expect(removeDay(df, c, 3), null);
    // day 0 taken back leaves the animals listed, with no day and no reading
    removeDay(df, c, 0);
    expect([df.rowCount, df.col('day')!.isNone(0), df.col('volume')!.isNone(1)].join(), '2,true,true');
  });

  test('two readings of one animal on one day are said as such; a removal reason takes an animal off', async () => {
    expect(twoReadings([{animal: 'A1', day: 3, volume: 180, weight: 20}, {animal: 'A1', day: 3, volume: 999,
      weight: 20}]), 'A1 has two volumes on day 3: keep one row per animal and day in the study table (correct or ' +
      'delete the other).');
    // a volume on one row and the weight on another is one reading of each
    expect(twoReadings([{animal: 'A1', day: 3, volume: 180, weight: null}, {animal: 'A1', day: 3, volume: null,
      weight: 20}]), null);
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['a1', 'a1', 'a2']),
      DG.Column.fromList(DG.TYPE.FLOAT, 'day', [0, 21, 21]), DG.Column.fromStrings('removal_reason', ['', 'tumor', '']),
      DG.Column.fromStrings('Mark reason', ['', '', 'typo'])]);
    const off = removedAnimals(df, df.col('animal')!, df.col('day')!);
    expect([...off].map(([a, x]) => `${a} ${x.day} ${x.why}`).join(), 'a1 21 tumor');
  });

  test('nothing is written until the day is saved, then every row at once, logged', async () => {
    const df = study();
    const c = columnsOf(df);
    const list = [animal('a1', 'Vehicle', 100), animal('a2', 'Drug', 120)];
    // an impossible value writes nothing at all
    let threw = false;
    try {
      writeDay(df, c, 3, list, {a1: {volume: '150'}, a2: {volume: '-4'}});
    } catch {
      threw = true;
    }
    expect(threw, true);
    expect(df.rowCount, 2);
    expect(logOf(df).length, 0);
    writeDay(df, c, 3, list, {a1: {volume: '150', weight: '20.5'}, a2: {missed: 'Dead'}});
    expect(df.rowCount, 3);
    expect(df.get('volume', 2), 150);
    expect(df.get('group', 2), 'Vehicle');
    expect(df.get('body_weight', 2), 20.5);
    expect(logOf(df).pop()!.what.startsWith('Day 3: 1 value entered; not measured: a2 (dead)'), true);
    expect(metaOf(df)!.notes.some((n) => n.animal === 'a2' && n.text === 'Not measured: dead'), true);
    // the day opened again starts from its values; a change edits its row, no new one
    expect(valuesOfDay(df, c, 3).a1.volume, '150');
    writeDay(df, c, 3, list, {a1: {volume: '160'}});
    expect(df.rowCount, 3);
    expect(df.get('volume', 2), 160);
  });

  test('a draft comes back for its study and day, and goes once saved', async () => {
    const key = `test-${Date.now()}`;
    saveDraft(key, 7, {a1: {volume: '150'}});
    expect(loadDraft(key, 7)?.values.a1.volume, '150');
    expect(loadDraft(key, 9), null);
    saveDraft(key, 7, null);
    expect(loadDraft(key, 7), null);
  });

  test('one reading corrected or flagged: logged, the animal untouched, read by Monitor and Analyze', async () => {
    const whole = await exampleStudy();
    const df = whole.clone(DG.BitSet.create(whole.rowCount, (i) => whole.get('day', i) <= 9));
    const inputs = {study: df, animal: df.col('animal')!, group: df.col('group')!, day: df.col('day')!,
      volume: df.col('volume')!, line: 'MC38', control: 'Vehicle', endpoint: 1500, correction: 'holm'};
    tagStudy({study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null});
    const n = {study: df, animal: 'animal', group: 'group', day: 'day', volume: 'volume', bodyWeight: null};
    const a = String(df.get('animal', 0));
    const row = lastReadingRow(n, a)!;
    const was = df.get('volume', row);
    const text = correctValues(df, row, {volume: 4321}, 'typo', `${a} on day ${df.get('day', row)}`);
    expect(text!.includes(`volume ${+was.toFixed(2)} → 4321 (typo)`), true);
    expect(logOf(df).pop()!.what, text);
    // the same value again is no change
    expect(correctValues(df, row, {volume: 4321}, 'typo', a), null);
    const frames = await monitorFrames(inputs as any);
    expect(frames.animals.get('last_volume', frames.animals.col('animal')!.toList().indexOf(a)), 4321);
    const read = analysisDeclarationOf(withoutExcluded(inputs) as any).data;
    expect(read.volume.some((v, i) => read.animal[i] === a && v === 4321), true);
    // flagged alone: only that row, the animal reads flagged, still in the analysis
    mark(df, [row], 'flagged', 'swelling', `${a} on day 9`);
    expect(animalMarks(n).get(a), 'flagged');
    expect(df.col('Flagged')!.toList().filter((x) => x).length, 1);
    expect(logOf(df).pop()!.what.startsWith(`Flagged ${a} on day 9: swelling`), true);
    expect(analysisDeclarationOf(withoutExcluded(inputs) as any).data.volume.length,
      analysisDeclarationOf(inputs as any).data.volume.length);
  });

  test('day 0 of animals listed before any measurement: written on their own rows', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', ['m1', 'm2']),
      DG.Column.fromStrings('group', ['', '']), DG.Column.float('day', 2), DG.Column.float('volume', 2),
      DG.Column.float('body_weight', 2)]);
    setMetaOn(df, emptyMeta());
    const first = (a: string): EntryAnimal => ({animal: a, group: '', lastDay: null, lastVolume: null,
      lastWeight: null, lower: null, upper: null, expected: '—'});
    writeDay(df, columnsOf(df), 0, [first('m1'), first('m2')], {m1: {volume: '90', weight: '20'}, m2: {volume: '110'}});
    expect(df.rowCount, 2);
    expect(df.get('day', 0), 0);
    expect(df.get('volume', 1), 110);
    expect(df.get('body_weight', 0), 20);
  });

  // QA 0.1.3 browser check: day 0 planned for another date than today was saved by Save day with no question
  test('day 0 saved for a date that is not today: asked at Save day, nothing written until answered', async () => {
    const button = (d: DG.Dialog, text: (s: string) => boolean) =>
      [...d.root.querySelectorAll<HTMLButtonElement>('button')].find((b) => text(b.textContent?.trim() ?? ''))!;
    for (const answer of ['keep', 'today', 'cancel'] as const) {
      // a followed plan as Monitor holds it (QA 0.1.3 r5): no rows yet, day 0's date only in the table's details;
      // the caller hands no date (Monitor's cached one was null), the dialog reads the table it writes to
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', []), DG.Column.fromStrings('group', []),
        DG.Column.float('day', 0), DG.Column.float('volume', 0), DG.Column.float('body_weight', 0)]);
      setMetaOn(df, {...emptyMeta(), start: '2099-01-05'});
      const start = () => metaOf(df)?.start ?? null;
      const first: EntryAnimal = {animal: 'm1', group: '', lastDay: null, lastVolume: null, lastWeight: null,
        lower: null, upper: null, expected: '—'};
      saveDraft('test-day0-date', 0, null);
      enterDayDialog({df, cols: columnsOf(df), study: 'test-day0-date', day: 0, animals: () => [first],
        warning: () => null, showGroup: false, saved: () => {}});
      const d = modalOf('enter-day')!;
      const field = d.root.querySelector<HTMLInputElement>('input[data-col="volume"]')!;
      field.value = '95';
      field.dispatchEvent(new Event('input'));
      button(d, (t) => t === 'Save day').click();
      await new Promise((r) => setTimeout(r, 50));
      const q = modalOf('day0-date');
      expect(!!q, true, 'Save day asks first');
      expect(df.rowCount, 0, 'nothing written before the answer');
      if (answer === 'keep')
        button(q!, (t) => t.startsWith('Keep')).click();
      else if (answer === 'today')
        button(q!, (t) => t === 'Day 0 is today').click();
      else
        button(q!, (t) => t.toLowerCase() === 'cancel').click();
      await new Promise((r) => setTimeout(r, 50));
      expect(df.rowCount, answer === 'cancel' ? 0 : 1, answer);
      expect(start() === '2099-01-05', answer !== 'today', `${answer}: day 0's date`);
      // Cancel closes the question only: the entry stays open with what was typed (QA 0.1.3 r5); an answer saves
      // the day and closes it
      expect(d.root.isConnected && modalOf('enter-day') === d, answer === 'cancel', `${answer}: the entry open`);
      if (answer === 'cancel')
        expect(d.root.querySelector<HTMLInputElement>('input[data-col="volume"]')!.value, '95');
      d.close();
      saveDraft('test-day0-date', 0, null);
    }
  });

  // QA 0.1.3 r5: a study with no animals in its table has their IDs only in the entry; closed with a draft, the
  // reopened entry rebuilt the default IDs and the values typed under M1…M9 were lost
  test('a day\'s draft keeps its animals\' IDs as typed, every value back under its own animal', async () => {
    const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('animal', []), DG.Column.fromStrings('group', []),
      DG.Column.float('day', 0), DG.Column.float('volume', 0), DG.Column.float('body_weight', 0)]);
    setMetaOn(df, emptyMeta());
    const animal = (a: string): EntryAnimal => ({animal: a, group: '', lastDay: null, lastVolume: null,
      lastWeight: null, lower: null, upper: null, expected: '—'});
    saveDraft('test-day0-ids', 0, null);
    enterDayDialog({df, cols: columnsOf(df), study: 'test-day0-ids', day: 0, animals: () => ['M1', 'M2', 'M10']
      .map(animal), warning: () => null, showGroup: false, saved: () => {}});
    const d = modalOf('enter-day')!;
    try {
      const field = d.root.querySelector<HTMLInputElement>('input[data-col="volume"]')!;
      field.value = '95';
      field.dispatchEvent(new Event('input'));
      // closed with Esc: asked, and the draft kept
      window.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true, cancelable: true}));
      await new Promise((r) => setTimeout(r, 50));
      modalOf('enter-day-close')?.close();
      const kept = loadDraft('test-day0-ids', 0);
      expect(JSON.stringify(kept?.animals), JSON.stringify(['M1', 'M2', 'M10']));
      expect(kept?.values.M1?.volume, '95');
    } finally {
      d.close();
      saveDraft('test-day0-ids', 0, null);
    }
  });
});
