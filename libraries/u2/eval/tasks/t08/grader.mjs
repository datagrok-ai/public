import {css, expect, flush, mount, need} from '../../lib/harness.mjs';

const CAMPAIGNS = [
  {code: 'KHL-12', name: 'PI3K quick wins', status: 'In Progress', created: new Date(2026, 9, 19)},
  {code: 'KHL-9', name: 'BRD4 follow-up', status: 'Submitted', created: new Date(2026, 8, 2)},
  {code: 'KHL-3', name: 'Kinase hits', status: 'In Progress', created: new Date(2026, 2, 30)},
];

async function build(mod) {
  const root = mount(need(mod, 'campaignsTable')(CAMPAIGNS.map((c) => ({...c}))));
  await flush();
  const items = root.querySelectorAll('[data-u2-item]');
  const sections = items.filter((el) => el.dataset.u2Item.startsWith('section:'));
  const rows = items.filter((el) => !el.dataset.u2Item.startsWith('section:') && el.tagName === 'TR');
  return {root, sections, rows};
}

const firstCell = (tr) => tr.querySelector('td, th') ?? tr;
const rgb = (c) => `rgb(${c.join(', ')})`;
const differ = (a, b) => a.some((x, i) => Math.abs(x - b[i]) >= 2);

export const checks = {
  'groups': async (mod) => {
    const {sections, rows} = await build(mod);
    const labels = sections.map((s) => s.textContent);
    expect(sections.length === 2, `${sections.length} group headers (${labels.join(' / ')}), not one per status`);
    expect(/In Progress.*2/.test(labels[0]) && /Submitted.*1/.test(labels[1]),
      `the group headers read ${labels.join(' / ')}, not each status with its count`);
    expect(rows.length === 3, `${rows.length} campaign rows under the groups`);
  },
  'hover-distinct': async (mod) => {
    const {sections, rows} = await build(mod);
    expect(sections.length > 0 && rows.length > 0, 'no group header and row to compare');
    const header = css.background(firstCell(sections[0]));
    const cells = rows[0].querySelectorAll('td, th');
    for (const [i, cell] of (cells.length ? cells : [rows[0]]).entries()) {
      const hovered = css.background(cell, {hover: cell});
      expect(differ(hovered, css.background(cell)),
        `cell ${i + 1} of the row under the pointer looks as at rest (${rgb(hovered)})`);
      expect(differ(hovered, header),
        `cell ${i + 1} of the hovered row has the group header's colour (${rgb(hovered)})`);
    }
  },
};
