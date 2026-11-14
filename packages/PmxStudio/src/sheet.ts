/* The randomization sheet: a page to print for the animal room, from the study table once its groups are set — the
   study, day 0 and the seed, what each group receives, then each animal with its group and its volume at
   randomization, a column to write the cage in, and the lines to sign. It shows the groups: printing it is logged. */
import * as DG from 'datagrok-api/dg';

import {roleColumns} from './study-file';
import {StudyMeta, dateOf, dateText} from './study-meta';

/** Text as HTML: a study's names (animals, groups, notes) come from files and are never markup. */
export const esc = (s: unknown) => String(s ?? '').replace(/[&<>"]/g, (c) =>
  ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'})[c]!);

/** Whether a measurement on `day` is closer to randomization than one on `was`: the last on or before day 0, else
 *  the first after it. */
export const closer = (day: number, was: number) => day <= 0 ? was > 0 || day > was : was > 0 && day < was;

/** Each animal once, by group then name, with its volume at randomization. */
export function sheetRows(df: DG.DataFrame): {animal: string, group: string, volume: number | null}[] {
  const c = roleColumns(df);
  const at = new Map<string, {group: string, day: number, volume: number | null}>();
  for (let i = 0; i < df.rowCount; i++) {
    const animal = String(c.animal.get(i));
    const day = c.day.get(i) as number;
    const was = at.get(animal);
    if (!was || closer(day, was.day)) {
      at.set(animal, {group: c.group.isNone(i) ? '' : String(c.group.get(i)), day,
        volume: c.volume.isNone(i) ? null : c.volume.get(i)});
    }
  }
  return [...at.entries()].map(([animal, x]) => ({animal, group: x.group, volume: x.volume}))
    .sort((a, b) => a.group.localeCompare(b.group) || a.animal.localeCompare(b.animal, undefined, {numeric: true}));
}

/** The sheet as a page of its own. */
export function sheetHtml(name: string, meta: StudyMeta | null, df: DG.DataFrame, seed: string | null): string {
  const rows = sheetRows(df);
  const day0 = dateOf(meta?.start ?? null, 0);
  const receives = (group: string) => {
    const g = meta?.groups.find((x) => x.group === group);
    return g ? [g.compound, g.dose === null ? '' : `${g.dose} ${g.unit}`, g.route, g.schedule].filter((x) => x)
      .join(' · ') : '';
  };
  const groups = [...new Set(rows.map((r) => r.group))];
  const sub = [meta?.id, day0 ? `Day 0: ${dateText(day0, true)}` : 'Day 0: ____________', seed ? `Seed ${seed}` : '',
    `${rows.length} animals`].filter((x) => x).map(esc).join(' · ');
  return `<!doctype html><html><head><meta charset="utf-8"><title>Randomization sheet · ${esc(name)}</title>
<style>
body{font:12px/1.4 system-ui,sans-serif;color:#222;margin:24px}
h1{font-size:18px;margin:0 0 4px}.sub{color:#555;margin-bottom:16px}
table{border-collapse:collapse;width:100%;margin-bottom:16px}
th,td{border:1px solid #999;padding:4px 6px;text-align:left}
th{background:#eee}td.n{text-align:right}
.sign{display:flex;gap:48px;margin-top:32px}.sign div{flex:1;border-top:1px solid #333;padding-top:4px}
@media print{body{margin:0}tr{break-inside:avoid}}
</style></head><body>
<h1>Randomization sheet · ${esc(name)}</h1>
<div class="sub">${sub}</div>
<table><tr><th>Group</th><th>Animals</th><th>Receives</th></tr>
${groups.map((g) => `<tr><td>${esc(g)}</td><td class="n">${rows.filter((r) => r.group === g).length}</td>` +
  `<td>${esc(receives(g))}</td></tr>`).join('\n')}</table>
<table><tr><th>Animal</th><th>Group</th><th>Volume at randomization (mm³)</th><th>Cage</th><th>Notes</th></tr>
${rows.map((r) => `<tr><td>${esc(r.animal)}</td><td>${esc(r.group)}</td><td class="n">` +
  `${r.volume === null ? '' : Math.round(r.volume)}</td><td></td><td></td></tr>`).join('\n')}</table>
<div class="sign"><div>Randomized by, date</div><div>Checked by, date</div></div>
</body></html>`;
}

/** Opens the sheet in a tab of its own, its print dialog up. */
export function printSheet(html: string): void {
  const w = window.open('', '_blank');
  if (!w)
    throw new Error('The browser blocked the sheet\'s tab: allow pop-ups for Datagrok and print again');
  w.document.write(html);
  // printed once its images are in. Listened for after the write: writing opens the document anew, and that drops
  // the listeners its window had
  w.addEventListener('load', () => {
    w.focus();
    w.print();
  });
  w.document.close();
}
