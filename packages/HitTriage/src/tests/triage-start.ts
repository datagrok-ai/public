/* eslint-disable max-len */
import {Scope} from '@datagrok-libraries/u2';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {newCampaignAccordeon} from '../app/accordeons/new-campaign-accordeon';
import {HitTriageDataSourceTag, PepTriageDataSourceTag, i18n} from '../app/consts';
import {HitTriageTemplate} from '../app/types';
import {expectSame, hasBalloon, pepTriageTemplate, Sandbox, triageTemplate} from './fixtures';

/** The new-campaign form of [template] with the inputs of the given captions typed in, START pressed: what START
 * handed on to create. */
async function startOver(template: HitTriageTemplate, tag: string, typed: {[caption: string]: string}): Promise<unknown[]> {
  const owner = new Scope();
  const started: unknown[] = [];
  try {
    const form = await newCampaignAccordeon(template, {}, tag, {owner, start: async (c) => {started.push(c);}});
    document.body.append(form);
    for (const [caption, value] of Object.entries(typed)) {
      const root = Array.from(form.querySelectorAll<HTMLElement>('.ui-input-root'))
        .find((r) => r.querySelector('label')?.textContent?.trim().toLowerCase() === caption.toLowerCase());
      const editor = root?.querySelector<HTMLInputElement>('input');
      expect(editor != null, true, `the ${caption} input`);
      editor!.value = value;
      editor!.dispatchEvent(new Event('input', {bubbles: true}));
      editor!.dispatchEvent(new Event('change', {bubbles: true}));
    }
    const start = Array.from(form.querySelectorAll('button')).find((b) => b.textContent?.trim().toLowerCase() === 'start')!;
    await waitUntil(() => !start.disabled, `START stayed disabled: ${start.title}`, 10000);
    start.click();
    await waitUntil(() => start.disabled, 'START never got busy', 5000).catch(() => {});
    await waitUntil(() => !start.disabled, 'START stayed busy', 30000);
    form.remove();
  } finally {
    owner.dispose();
  }
  return started;
}

/** Each caption's right edge in the new-campaign section of a landing, by its text; START's and the form's right edges. */
function geometry(root: HTMLElement): {labels: {[text: string]: number}, start: number, form: number} {
  const picker = root.querySelector<HTMLElement>('[data-u2-name="template-picker"]')!;
  const labels: {[text: string]: number} = {};
  for (const label of Array.from(picker.querySelectorAll<HTMLElement>('label, .u2-input-label'))) {
    if (label.textContent?.trim() && label.getClientRects().length > 0)
      labels[label.textContent.trim()] = Math.round(label.getBoundingClientRect().right);
  }
  const start = Array.from(picker.querySelectorAll('button')).find((b) => b.textContent?.trim().toLowerCase() === 'start');
  const form = picker.querySelector<HTMLElement>('.ui-form');
  return {labels, start: Math.round(start?.getBoundingClientRect().right ?? 0), form: Math.round(form?.getBoundingClientRect().right ?? 0)};
}

category('New campaign: Hit Triage, PepTriage', () => {
  test('Hit Triage: Template, Campaign Name and the data source share one label column, START at the form\'s end', async () => {
    const sb = await Sandbox.create('htform');
    try {
      await sb.writeTemplate(triageTemplate('ZZB3 HT', 'ZZB3'));
      const app = await sb.hitTriage();
      const root = app.infoView.root;
      let g = geometry(root);
      const near = (a?: number, b?: number) => a !== undefined && b !== undefined && Math.abs(a - b) <= 1;
      await waitUntil(() => (g = geometry(root)).start > 0 && near(g.labels['Template'], g.labels['Campaign Name']),
        () => `Template does not line up with Campaign Name: ${JSON.stringify(g)}`, 10000);
      expect(near(g.labels[i18n.dataSourceFunction], g.labels['Campaign Name']), true, `the data source caption: ${JSON.stringify(g)}`);
      expect(Math.abs(g.start - g.form) <= 2, true, `START ends where the form does: ${JSON.stringify(g)}`);
    } finally {
      await sb.dispose();
    }
  }, {timeout: 60000});

  test('PepTriage: a data source that answers no rows creates no campaign and says so', async () => {
    const started = await startOver(pepTriageTemplate('ZZB3 PT', 'ZZB3'), PepTriageDataSourceTag, {Chemist: 'Ann'});
    expectSame(started, [], 'the campaigns START created');
    expect(hasBalloon('Demo Peptide Sequences returned no rows, so no campaign was created'), true, 'the warning');
  }, {timeout: 60000});

  test('Hit Triage: a data source that answers no rows creates no campaign and says so', async () => {
    const started = await startOver(triageTemplate('ZZB3 HT', 'ZZB3'), HitTriageDataSourceTag,
      {'Number Of Molecules': '0', 'Chemist': 'Ann'});
    expectSame(started, [], 'the campaigns START created');
    expect(hasBalloon('Demo Molecules variable returned no rows, so no campaign was created'), true, 'the warning');
  }, {timeout: 60000});
});
