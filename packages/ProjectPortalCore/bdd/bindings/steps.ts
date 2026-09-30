/* The steps only this package can define: how to reach the app, a scratch program to take
   through a gate, and reading a board's stage columns in order. */
import type {Page} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

export const APP_PATH = '/apps/ProjectPortalCore';

export const openApp = Given('user opens the Project Portal app', async (page: Page) => {
  await page.goto(APP_PATH, {waitUntil: 'domcontentloaded'});
  await page.locator('[data-u2-name="portfolio-list"] .pp-portfolio').first().waitFor({timeout: 60000});
}, {tier: 'ui', enters: 'Project Portal app',
  description: 'needs this package published and ProjectPortalCore:loadPortfolios + syncAllPortfolios run on the stand'});

/** Removes every live program with [code] in [portfolio]; the sweep at the start and the cleanup at the end. */
async function removeScratch(page: Page, code: string, portfolio: string): Promise<number> {
  return page.evaluate(async ([c, p]) => {
    const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
    const programs = grok.dapi.domains.table('portal.program');
    const rows = await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`});
    for (const r of rows)
      await programs.delete(r.id);
    return (await programs.query({filter: `portfolio_id = "${pf.id}" and code = "${c}"`})).length;
  }, [code, portfolio]);
}

export const scratchProject = Given('a scratch project {string} in the {string} portfolio',
  async (page: Page, code: string, portfolio: string) => {
    silent(page);   // set-up a person does not do: kept out of guide videos
    await removeScratch(page, code, portfolio);
    atFeatureEnd(page, async () => {
      const left = await removeScratch(page, code, portfolio);
      if (left !== 0)
        throw new Error(`${left} scratch program(s) ${code} survived the cleanup`);
    });
    await page.evaluate(async ([c, p]) => {
      const pf = await grok.dapi.domains.table('portal.portfolio').getByKey({code: p});
      await grok.dapi.domains.table('portal.program').insert({portfolio_id: pf.id, code: c, name: c,
        stage_model_id: pf.stage_model_id, source: 'manual'});
    }, [code, portfolio]);
    // the board was drawn before the program existed: redraw it
    await page.evaluate(async () => {
      const view = grok.shell.v;
      await view?.portal?.load?.();
    });
  }, {tier: 'api', description: 'a manual program, swept first and deleted at the end of the feature'});

export const boardStages = Then('the board should show the stages {string}', async (page: Page, codes: string) => {
  const expected = codes.split(',').map((s) => s.trim());
  await expect.poll(() => page.locator('[data-u2-name="stage-board"] .pp-column-title').allTextContents()
    .then((titles) => titles.map((t) => t.split(' · ')[0]).filter((c) => c !== 'Unstaged')),
  {message: `stage columns in order`}).toEqual(expected);
}, {description: 'the stage codes of the selected portfolio\'s model, left to right, unstaged column excluded'});

export const sidePanelsHidden = Given('the side panels are hidden', async (page: Page) => {
  silent(page);   // framing for a guide, not something the viewer does
  const before = await page.evaluate(() => {
    const w = grok.shell.windows;
    const was = {browse: w.showBrowse, context: w.showContextPanel, help: w.showHelp};
    w.showBrowse = false;
    w.showContextPanel = false;
    w.showHelp = false;
    // a dev stand announces its debug-mode packages; not part of the story
    document.querySelectorAll('.d4-balloon').forEach((b) => b.remove());
    return was;
  });
  atFeatureEnd(page, async () => {
    await page.evaluate((was) => {
      grok.shell.windows.showBrowse = was.browse;
      grok.shell.windows.showContextPanel = was.context;
      grok.shell.windows.showHelp = was.help;
    }, before);
  });
}, {tier: 'ui', description: 'hides the browse tree, context panel and help for the rest of the feature'});
