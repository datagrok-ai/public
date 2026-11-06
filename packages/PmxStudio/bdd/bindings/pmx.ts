/* PMX Studio's own names: the form's Run (the Context Panel shows the app's Run beside it), the wait for a plan or an
   analysis, read from the view's loader, and the study's figures, remembered and compared as a change of the form
   should move them (the numbers themselves are the package tests'). */
import {Page, test} from '@playwright/test';
import {element, Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, openTopMenu, pollMs, silent} from '@datagrok-libraries/bdd/runtime';

element('form run button', {selector: '.pmx-declaration [data-u2="button"]:has-text("Run")'});
element('advanced header', {selector: '.pmx-declaration .u2-form-category:has-text("Advanced")'});
element('stale note', {selector: '.pmx-stale'});
element('study error', {selector: '.pmx-plan-error, .pmx-analyze-error'});
element('calendar status', {selector: '.pmx-calendar-status'});
element('calendar summary', {selector: '.pmx-calendar-summary'});
element('study calendar', {selector: '.pmx-calendar'});
element('analyze start', {selector: '.pmx-analyze-start'});
element('animal list', {selector: '.pmx-animal-list'});
element('twin effect line', {selector: '.pmx-twin-side'});
element('results pane', {selector: '.pmx-results'});
element('tc choice', {selector: '.pmx-arm .u2-radio-item:has-text("T/C")'});
element('linear axis choice', {selector: '.u2-radio-item:visible:has-text("Linear")'});
element('log axis choice', {selector: '.u2-radio-item:visible:has-text("Log")'});
element('visible volume axis', {selector: '[data-u2="radio-input"]:visible:has(.u2-radio-item:has-text("Linear"))'});
element('tuesday header', {selector: '.pmx-calendar-weekday:has-text("Tue")'});
element('thursday header', {selector: '.pmx-calendar-weekday:has-text("Thu")'});
element('randomization day picker', {selector: '.pmx-calendar-summary [data-u2="choice-input"]'});
element('first copy button', {selector: '[data-u2="button"]:visible:has-text("Copy") >> nth=0'});
// the redesign: the Overview's cards, the stage stepper, the layout editor's tiles and the Analyze charts
element('plan card', {selector: '.pmx-card[data-stage="plan"]'});
element('analyze step', {selector: '.pmx-head-steps button:has-text("Analyze")'});
element('growth rates tile', {selector: '.pmx-layout-tile[data-panel="rates"]'});
element('growth rates corner', {selector: '.pmx-layout-tile[data-panel="rates"] .pmx-layout-grip'});
element('survival tile', {selector: '.pmx-layout-tile[data-panel="km"]'});
element('waterfall chart', {selector: '.pmx-waterfall-chart'});
element('first animal row', {selector: '.pmx-forest-row >> nth=0'});
// the row, not the input: the drag step looks for the range input inside the element
element('km slider', {selector: '.pmx-km-row'});
element('km chart', {selector: '.pmx-km-chart'});
element('results tab', {selector: '.pmx-tab:has-text("Results")'});
element('each animal tab', {selector: '.pmx-tab:has-text("Each animal")'});

const figure = (page: Page, label: string) => page.locator('[data-u2="stat-card"]').filter({visible: true})
  .filter({hasText: label}).first();
/** A figure's value: the card's text without its label. */
const reading = async (page: Page, label: string): Promise<string> =>
  (await figure(page, label).innerText()).replace(label, '').trim();
const amount = (text: string) => Number(/-?\d+(\.\d+)?/.exec(text)?.[0] ?? NaN);

const REMEMBERED = new WeakMap<Page, Map<string, string>>();

export const finishedRunning = Then('the study should have finished running', async (page: Page) => {
  silent(page);
  await expect.poll(() => page.locator('.pmx-busy').filter({visible: true}).count(),
    {timeout: pollMs(120000), message: 'the view\'s loader shows'}).toBe(0);
}, {description: 'the view\'s loader is gone: the plan or the analysis is on the page; silent in a guide'});

export const effectsArrived = Then('the smallest detectable effects should have arrived', async (page: Page) => {
  const cell = page.locator('tbody tr').filter({visible: true, hasText: 'Treatment 1'}).locator('td').last();
  await expect.poll(() => cell.innerText(), {timeout: pollMs(120000), message: 'the arm\'s smallest detectable effect'})
    .toMatch(/^(TGI|T\/C) \d+%$|^out of reach$/);
}, {description: 'with the animals per arm set, the groups table\'s last column is filled in after the plan'});

export const rememberFigure = When('user remembers the {string} figure', async (page: Page, label: string) => {
  const map = REMEMBERED.get(page) ?? new Map<string, string>();
  map.set(label, await reading(page, label));
  REMEMBERED.set(page, map);
}, {tier: 'api', description: 'the figure as the card shows it, for a later comparison'});

const remembered = (page: Page, label: string): string => {
  const was = REMEMBERED.get(page)?.get(label);
  if (was === undefined)
    throw new Error(`the "${label}" figure was not remembered`);
  return was;
};

export const figureHigher = Then('the {string} figure should be higher than remembered',
  async (page: Page, label: string) => {
    const was = amount(remembered(page, label));
    await expect.poll(async () => amount(await reading(page, label)),
      {message: `the "${label}" figure, ${was} before`}).toBeGreaterThan(was);
  });

export const figureLower = Then('the {string} figure should be lower than remembered',
  async (page: Page, label: string) => {
    const was = amount(remembered(page, label));
    await expect.poll(async () => amount(await reading(page, label)),
      {message: `the "${label}" figure, ${was} before`}).toBeLessThan(was);
  });

export const figureSame = Then('the {string} figure should be as remembered', async (page: Page, label: string) => {
  const was = remembered(page, label);
  await expect.poll(() => reading(page, label), {message: `the "${label}" figure`}).toBe(was);
});

export const figureReads = Then('the {string} figure should read {string}',
  async (page: Page, label: string, text: string) => {
    await expect.poll(() => reading(page, label), {message: `the "${label}" figure`}).toBe(text);
  });

export const noSavedPlan = Given('no plan named {string} is saved', async (page: Page, name: string) => {
  const remove = (n: string) => page.evaluate(async (plan) => {
    const g = (window as any).grok;
    const path = `${g.shell.user.project.name}:Home/PMX Studio/plans/${plan}.pmxplan`;
    if (await g.dapi.files.exists(path))
      await g.dapi.files.delete(path);
    return g.dapi.files.exists(path);
  }, n);
  expect(await remove(name), 'the plan is gone before the feature').toBe(false);
  atFeatureEnd(page, async () => expect(await remove(name), 'the plan is gone after the feature').toBe(false));
}, {tier: 'api', description: 'removed now and at feature end, and read back'});

export const planIsSaved = Then('a plan named {string} should be saved', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(async (plan) => {
    const g = (window as any).grok;
    return g.dapi.files.exists(`${g.shell.user.project.name}:Home/PMX Studio/plans/${plan}.pmxplan`);
  }, name), {message: `the plan "${name}" in the user's plans folder`}).toBe(true);
}, {description: 'read back from the server'});

export const pickSavedPlan = When('user picks the saved plan {string}', async (page: Page, name: string) => {
  await page.locator('.pmx-declaration [data-u2$="-input"]:visible').filter({hasText: 'Plan'})
    .locator('.fa-list, [name*="list" i]').first().click();
  await page.locator('.d4-menu-item').filter({hasText: `${name}.pmxplan`}).first().click();
}, {tier: 'ui', description: 'through the list icon of the Plan field in Analyze'});

export const rememberTwin = When('user remembers the twin line', async (page: Page) => {
  const map = REMEMBERED.get(page) ?? new Map<string, string>();
  map.set('twin', await page.locator('.pmx-twin-side').innerText());
  REMEMBERED.set(page, map);
}, {tier: 'api'});

export const pickLastAnimal = When('user clicks on the last animal row', async (page: Page) => {
  await page.locator('.pmx-forest-row').last().click();
}, {tier: 'ui'});

export const twinDiffers = Then('the twin line should differ from remembered', async (page: Page) => {
  await expect.poll(() => page.locator('.pmx-twin-side').innerText(), {message: 'the animal drawn'})
    .not.toBe(remembered(page, 'twin'));
});

// the teaser guide: the stage steps, Monitor's attention list, the Overview's studies, and what filming needs
element('overview step', {selector: '.pmx-head-steps button:has-text("Overview")'});
element('monitor step', {selector: '.pmx-head-steps button:has-text("Monitor")'});
element('first attention row', {selector: '.pmx-monitor-attention .u2-table-row >> nth=0'});
element('your studies', {selector: '.pmx-library'});
// the demo MC38 study at its three moments, as the Overview lists it
element('demo plan row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response (its plan)")'});
element('demo followed row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response"):not(:has-text("(its plan)")):not(:has-text("(all its data)"))'});
element('demo finished row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response (all its data)")'});

export const presenting = Given('the shell is in presentation mode without tooltips', async (page: Page) => {
  silent(page);
  // a film shows the app, not the platform's panels, and no hint a pointer left over a control would raise
  await page.addStyleTag({content: '.d4-tooltip, [role="tooltip"] { display: none !important; }'});
  await page.evaluate(() => {
    const strip = (root: ParentNode) => root.querySelectorAll('[title]').forEach((e) => e.removeAttribute('title'));
    strip(document);
    new MutationObserver((records) => records.forEach((r) => {
      if (r.type === 'attributes')
        (r.target as Element).removeAttribute('title');
      else
        r.addedNodes.forEach((n) => n instanceof Element && (n.removeAttribute('title'), strip(n)));
    })).observe(document.body, {subtree: true, childList: true, attributes: true, attributeFilter: ['title']});
    (window as any).grok.shell.windows.presentationMode = true;
  });
}, {tier: 'api', description: 'for a guide that is filmed: presentation mode on, tooltips and native titles hidden'});

export const pointAtMenu = When('user points at {string} in the top menu', async (page: Page, path: string) => {
  await openTopMenu(page, path, false);
  await page.waitForTimeout(1500);
  await page.keyboard.press('Escape');
  await page.mouse.move(960, 600);
}, {tier: 'ui', description: 'the menu walked to the item and left open a moment, nothing run'});

export const stageOpen = Then('the {word} stage should be open', async (page: Page, stage: string) => {
  silent(page);
  await page.locator('.pmx-head-step.on').filter({visible: true, hasText: new RegExp(`^${stage}$`, 'i')}).first()
    .waitFor({timeout: pollMs(120000)});
  await expect.poll(() => page.locator('.pmx-busy').filter({visible: true}).count(),
    {timeout: pollMs(120000), message: 'the view\'s loader shows'}).toBe(0);
}, {description: 'the stage\'s view is current and its loader gone (a demo study is simulated as it opens); silent in a guide'});

export const guideBudget = Given('the guide may take up to {int} minutes', async (page: Page, minutes: number) => {
  test.info().setTimeout(Math.max(test.info().timeout, minutes * 60000));
}, {tier: 'api', description: 'a filmed walk through every stage outlasts a test\'s budget'});

export const studiesListed = Then('the studies should be listed', async (page: Page) => {
  silent(page);
  await expect.poll(() => page.locator('.pmx-library .u2-table-row').filter({visible: true}).count(),
    {timeout: pollMs(120000), message: 'the Overview\'s studies'}).toBeGreaterThan(3);
}, {description: 'the Overview lists the user\'s and the demo studies (each demo is simulated first); silent in a guide'});
