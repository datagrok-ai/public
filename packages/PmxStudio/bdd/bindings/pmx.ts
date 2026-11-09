/* PMX Studio's own names: the form's Run (the Context Panel shows the app's Run beside it), the wait for a plan or an
   analysis, read from the view's loader, and the study's figures, remembered and compared as a change of the form
   should move them (the numbers themselves are the package tests'). */
import {Page, test} from '@playwright/test';
import {element, Given, lookupElement, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, gestures, pollMs, silent, unlit} from '@datagrok-libraries/bdd/runtime';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';

element('form run button', {selector: '.pmx-declaration [data-u2="button"]:has-text("Run")'});
element('advanced header', {selector: '.pmx-declaration .u2-form-category:has-text("Advanced")'});
element('stale note', {selector: '.pmx-stale'});
element('study error', {selector: '.pmx-stage-error'});
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
element('analyze header step', {selector: '.pmx-head-steps button:has-text("Analyze")'});
element('first attention row', {selector: '.pmx-monitor-attention .u2-table-row >> nth=0'});
element('first sheet row', {selector: '.pmx-monitor-sheet .u2-table-row >> nth=0'});
element('your studies', {selector: '.pmx-library'});
// what the teaser points out without touching: Plan's form and figures, Analyze's twin
element('plan figures', {selector: '.pmx-plan-cards'});
// a section by its title: the innermost one is the last of those that hold it
element('tumor model section', {selector: '.u2-section:has(.u2-section-title:text-is("Tumor model"))'});
element('treatment arms section', {selector: '.u2-section:has(.u2-section-title:text-is("Treatment arms"))'});
element('experimental groups section', {selector: '.u2-section:has(.u2-section-title:text-is("Experimental groups"))'});
element('twin section', {selector: '.u2-section:has(.u2-section-title:text-is("Each animal against its untreated twin"))'});
element('tolerability tab', {selector: '.pmx-tab:has-text("Tolerability")'});
element('third animal row', {selector: '.pmx-forest-row >> nth=2'});
// the demo MC38 study at its three moments, as the Overview lists it
element('demo plan row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response (its plan)")'});
element('demo followed row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response"):not(:has-text("(its plan)")):not(:has-text("(all its data)"))'});
element('demo finished row', {selector: '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response (all its data)")'});

/** For a guide that is filmed: no hint a pointer left over a control would raise (tooltips, native titles), nor the
 *  hover tint of a row nothing happens on (Monitor's lists), which reads as a row picked out. */
async function noHints(page: Page): Promise<void> {
  await page.addStyleTag({content: '.d4-tooltip, [role="tooltip"] { display: none !important; } ' +
    ':is(.pmx-monitor-sheet, .pmx-monitor-attention) .u2-table-row:hover { background: var(--u2-row-base-bg) !important; } ' +
    ':is(.pmx-monitor-sheet, .pmx-monitor-attention) .u2-table-row:hover > :first-child { box-shadow: none !important; }'});
  await page.evaluate(() => {
    const strip = (root: ParentNode) => root.querySelectorAll('[title]').forEach((e) => e.removeAttribute('title'));
    strip(document);
    new MutationObserver((records) => records.forEach((r) => {
      if (r.type === 'attributes')
        (r.target as Element).removeAttribute('title');
      else
        r.addedNodes.forEach((n) => n instanceof Element && (n.removeAttribute('title'), strip(n)));
    })).observe(document.body, {subtree: true, childList: true, attributes: true, attributeFilter: ['title']});
  });
}

export const noTooltips = Given('the shell shows no tooltips', async (page: Page) => {
  silent(page);
  await noHints(page);
}, {tier: 'api', description: 'for a guide that is filmed: tooltips and native titles hidden; silent in a guide'});

export const presenting = Given('the shell is in presentation mode without tooltips', async (page: Page) => {
  silent(page);
  // a film shows the app, not the platform's panels
  await noHints(page);
  await page.evaluate(() => {
    (window as any).grok.shell.windows.presentationMode = true;
  });
}, {tier: 'api', description: 'for a guide that is filmed: presentation mode on, tooltips and native titles hidden'});

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

/** What a guide films as the step's last picture is the page once PMX Studio has finished: a view a click opens (added
 *  a task later; a demo study is simulated first), and the loader of what it computes. */
async function settled(page: Page, viewChanges: boolean, before: string): Promise<void> {
  if (viewChanges)
    await page.waitForFunction((b) => (window as any).grok.shell.v?.name !== b, before, {timeout: pollMs(120000)});
  await page.waitForTimeout(400);
  // the view's loader, and the inline ones of what follows later (the twins in Analyze) — those on screen: a pane below
  // the fold (Plan's sweeps) is not in the picture
  await expect.poll(() => page.evaluate(() => Array.from(document.querySelectorAll('.pmx-busy, .u2-loader-spinner'))
    .filter((e) => {
      const r = e.getBoundingClientRect();
      return r.width > 0 && r.height > 0 && r.bottom > 0 && r.top < window.innerHeight;
    }).length), {timeout: pollMs(120000), message: 'a loader shows on screen'}).toBe(0);
}
const currentView = (page: Page) => page.evaluate(() => String((window as any).grok.shell.v?.name ?? ''));

export const clickOpensStage = When('user clicks on {element} and the stage opens', async (page: Page, target: ElementRef) => {
  const before = await currentView(page);
  await clearPointed(page);
  await gestures.click(page, target);
  await settled(page, true, before);
}, {tier: 'ui', description: 'the click, then the view it opens current and its loader gone'});

export const clickSettles = When('user clicks on {element} and PMX Studio settles', async (page: Page, target: ElementRef) => {
  await clearPointed(page);
  await gestures.click(page, target);
  await settled(page, false, '');
}, {tier: 'ui', description: "the click, then the view's loader gone"});

export const enterSettles = When('user enters {string} into {element} and PMX Studio settles',
  async (page: Page, text: string, target: ElementRef) => {
    await gestures.typeInto(page, target, text, true);
    await settled(page, false, '');
  }, {tier: 'ui', description: 'the value typed and committed, then the recompute it starts finished'});

export const openOverview = Given('user opens PMX Studio on its Overview', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => {
    const g = (window as any).grok;
    g.shell.addView(await g.functions.call('PmxStudio:pmxStudioApp', {}));
  });
  await expect.poll(() => page.locator('.pmx-library .u2-table-row').filter({visible: true}).count(),
    {timeout: pollMs(120000), message: "the Overview's studies"}).toBeGreaterThan(3);
  await page.mouse.move(12, 1070);
}, {tier: 'api', description: 'the app opened and its studies listed (each demo is simulated first); silent in a guide'});

export const studiesShown = Then("the lab's studies should be on the Overview", async (page: Page) => {
  await expect(page.locator('.pmx-library .u2-table-row').filter({visible: true}).first()).toBeVisible();
}, {description: 'the step a guide opens on: the Overview with its studies'});

/** Opens the demo MC38 study in Monitor (followed) and Analyze (all its data) from the Overview's list, each loaded,
 *  and comes back to the Overview: the header's stages then go to a study already there. */
export const demoOpen = Given('the demo study is open in Monitor and Analyze', async (page: Page) => {
  silent(page);
  const DEMO = '.pmx-library .u2-table-row:has-text("MC38 · Cpd-A dose response")';
  for (const [row, stage] of [[`${DEMO}:not(:has-text("(its plan)")):not(:has-text("(all its data)"))`, 'Monitor'],
    [`${DEMO}:has-text("(all its data)")`, 'Analyze']]) {
    await page.locator(row)
      .filter({visible: true}).first().click();
    await page.locator('.pmx-head-step.on').filter({visible: true, hasText: new RegExp(`^${stage}$`, 'i')}).first()
      .waitFor({timeout: pollMs(120000)});
    await settled(page, false, '');
    await page.locator('.pmx-head-steps button:has-text("Overview")').filter({visible: true}).first().click();
    await page.locator('.pmx-library .u2-table-row').filter({visible: true}).first().waitFor({timeout: pollMs(30000)});
  }
  await page.mouse.move(12, 1070);
}, {tier: 'api', description: 'the demo study opened in Monitor and Analyze and loaded, the Overview current again; silent in a guide'});

/** The element a selector registered for `target` names, the innermost one on screen (a section holds the sections
 *  inside it). Found by its selector, not through the runtime's lookup, so a guide does not light it as the step's
 *  target. */
function plainly(page: Page, target: ElementRef) {
  const selector = lookupElement(target.phrase)?.selector;
  if (!selector)
    throw new Error(`${target.phrase} has no selector`);
  return page.locator(selector).filter({visible: true}).last();
}

/** The yellow call-out on the page, if any; the next action clears it. */
const clearPointed = (page: Page) => page.evaluate(() => (window as any).bddPointed?.());

/** A yellow frame over the element until the next action: laid over the whole page, so no pane clips it and no field
 *  covers it, and following the element as it moves. Only looked at: the pointer stays where it is. */
export const pointOut = When('user points out {element}', async (page: Page, target: ElementRef) => {
  const el = plainly(page, target);
  await el.waitFor({timeout: pollMs(15000)});
  await clearPointed(page);
  await el.evaluate((e) => {
    const frame = document.createElement('div');
    Object.assign(frame.style, {position: 'fixed', border: '3px solid #F5C400', borderRadius: '8px',
      pointerEvents: 'none', zIndex: '100000', boxSizing: 'border-box'});
    document.body.append(frame);
    let on = true;
    const follow = () => {
      if (!on)
        return;
      const r = e.getBoundingClientRect();
      const shown = e.isConnected && r.width > 0 && r.height > 0;
      Object.assign(frame.style, {display: shown ? '' : 'none', left: `${r.left - 5}px`, top: `${r.top - 5}px`,
        width: `${r.width + 10}px`, height: `${r.height + 10}px`});
      requestAnimationFrame(follow);
    };
    follow();
    (window as any).bddPointed = () => {
      on = false;
      frame.remove();
      (window as any).bddPointed = undefined;
    };
  });
}, {tier: 'ui', description: "a yellow frame over it until the next action (a guide's picture shows it); nothing pressed"});

/** The pane scrolled by the wheel until the element is in view (in the middle, or its top near the pane's when it is
 *  taller): a guide films the scroll. The yellow frames go first, and the step is not lit. */
export const scrollTo = When('user scrolls the page to {element}', async (page: Page, target: ElementRef) => {
  unlit(page);
  await clearPointed(page);
  const el = plainly(page, target);
  await el.waitFor({timeout: pollMs(15000)});
  const way = await el.evaluate((e) => {
    let pane = e.parentElement;
    while (pane && !(pane.scrollHeight > pane.clientHeight && /auto|scroll/.test(getComputedStyle(pane).overflowY)))
      pane = pane.parentElement;
    if (!pane)
      return null;
    const [p, r] = [pane.getBoundingClientRect(), e.getBoundingClientRect()];
    const want = r.height < p.height * 0.8 ? (r.top + r.height / 2) - (p.top + p.height / 2) : r.top - p.top - 16;
    const dy = Math.max(-pane.scrollTop, Math.min(pane.scrollHeight - pane.clientHeight - pane.scrollTop, want));
    // the pane's own margin, not a chart in it that could take the wheel for itself
    return {x: p.left + 6, y: p.top + p.height / 2, dy};
  });
  if (!way)
    throw new Error(`${target.phrase} is in no scrolling pane`);
  if (Math.abs(way.dy) < 4)
    return;
  // the wheel turns where the pointer is: over the pane, near its left edge
  await page.mouse.move(way.x, way.y, {steps: 4});
  await page.mouse.wheel(0, way.dy);
  await page.waitForTimeout(400);
  const box = await el.boundingBox();
  const view = page.viewportSize();
  if (!box || !view || box.y > view.height - 40 || box.y + box.height < 40)
    throw new Error(`the wheel left ${target.phrase} out of view`);
}, {tier: 'ui', description: 'the pane wheeled until the element is in view; its yellow frames cleared first, not lit'});

/** A layout someone saved (a resized or moved panel) is not the one the stages define. */
export const defaultLayouts = Given('the stages show their default layouts', async (page: Page) => {
  silent(page);
  await page.evaluate(async () => {
    const settings = (window as any).grok.userSettings;
    for (const key of Object.keys(settings.get('PmxStudio') ?? {}).filter((k) => k.startsWith('layout-')))
      settings.delete('PmxStudio', key);
    await settings.flush();
  });
}, {tier: 'api', description: "the user's saved panel layouts removed; silent in a guide"});

export const dragSettles = When('user drags the slider of {element} to {float} and PMX Studio settles',
  async (page: Page, target: ElementRef, value: number) => {
    // what the drag changes is pointed out already (a yellow call-out): the step's pictures are not lit over it
    unlit(page);
    await gestures.dragSlider(page, target, value);
    await settled(page, false, '');
  }, {tier: 'ui', description: 'the slider dragged to the value, then the recompute it starts finished'});
