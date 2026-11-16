// Reactions in Crux mode, for the demo clips of features/demos/reactions-crux.feature: the showcase table opened from
// the project's fixtures (nothing written to the server), the Reaction pane's atom mapping waited for, the engine in
// effect read after a switch in Chem's settings (the Reaction pane is widgets/reaction-panel.ts; its template notes
// are an element of elements.ts). What Crux draws is Chem's package tests ('Reactions: Crux').
import {readFileSync} from 'node:fs';
import {basename, resolve} from 'node:path';
import {Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {el, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

export const openProjectCsv = Given('user opens the {string} file of the project', async (page: Page, file: string) => {
  const csv = readFileSync(resolve(process.env.BDD_ROOT ?? process.cwd(), file), 'utf8');
  const name = basename(file).replace(/\.[^.]+$/, '');
  await page.evaluate(async ([text, n]) => {
    const df = DG.DataFrame.fromCsv(text);
    df.name = n;
    // typed before the view opens, as a file is: the renderers are the semantic types' from the first paint
    await grok.data.detectSemanticTypes(df);
    const tv = grok.shell.addTableView(df);
    // the view makes row 0 of its first column current a second after the grid is made; done now, so that no viewer
    // repaints its marker mid-clip (as the library's dataset step does)
    const first = tv.grid?.columns.byIndex(1)?.column ?? df.columns.byIndex(0);
    df.currentCell = df.cell(0, first.name);
  }, [csv, name] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid &&
    !grok.shell.tv.grid.isRenderPending, name), {message: `a table view of "${name}", drawn`, timeout: pollMs(60000)}).toBe(true);
}, {tier: 'api', description: 'a CSV of the bdd project (`fixtures/…`) opened in a table view of its own, named as the file, its semantic types detected first; nothing is written to the stand'});

export const routeMapped = Then('the Reaction pane should have mapped {int} steps', async (page: Page, steps: number) => {
  const pane = await locate(page, el('"Reaction" section'));
  // RXNMapper maps a step at a time on a Crux worker, the first after its weights are read (seconds each)
  await expect(pane.first(), `step ${steps} of the route mapped`).toContainText(`Step ${steps}: mapped`,
    {timeout: pollMs(180000)});
  await expect(pane.first(), 'every step of the route mapped').not.toContainText('not mapped');
  await expect(pane.first(), 'no step given up').not.toContainText('took longer');
}, {description: 'the Reaction pane\'s Mapping status says every step of the route mapped (the panel maps by itself since 2026-11-14; polled for as long as RXNMapper takes); none "not mapped", none given up'});

export const scrollReactionDialog = When('user scrolls the Reaction dialog to its last step', async (page: Page) => {
  const dialog = (await locate(page, el('"Reaction" dialog'))).first();
  const first = dialog.locator('canvas').first();
  await expect(first, 'the full-screen picture').toBeVisible();
  const box = (await first.boundingBox())!;
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  // the full-screen view scrolls its pictures (the route, then its steps) in a box of its own: the wheel turns over it,
  // a notch at a time, until the last step is in view
  const left = () => first.evaluate((c) => {
    const host = c.parentElement!;
    return host.scrollHeight - host.clientHeight - host.scrollTop;
  });
  for (let i = 0; i < 8; i++) {
    const before = await left();
    if (before <= 1)
      break;
    await page.mouse.wheel(0, 300);
    await expect.poll(left, {message: 'the full-screen view scrolled by a notch'}).toBeLessThan(before);
  }
  expect(await left(), 'the full-screen view scrolled to its last step').toBeLessThanOrEqual(1);
}, {tier: 'ui', description: 'the wheel turned over the full-screen Reaction dialog a notch at a time until its last step picture is in view (a route\'s steps are under its picture)'});

export const pickOwnSettingsGroup = When('user picks their own group from the package settings groups', async (page: Page) => {
  const own = await page.evaluate(() => String(grok.shell.user.group?.friendlyName ?? grok.shell.user.login));
  await gestures.click(page, el(`${JSON.stringify(own)} package settings group`));
  await expect(await locate(page, el('package settings group button')), 'the group the Settings form edits')
    .toContainText(own);
}, {tier: 'ui', description: 'the signed-in account\'s own group, picked in the list the package settings group button opened; the form then edits its settings'});

export const engineInEffect = Then('the Chemical Engine should be {string}', async (page: Page, engine: string) => {
  await expect.poll(() => page.evaluate(() => grok.functions.call('Chem:chemicalEngine')),
    {message: 'the Chemical Engine this page runs', timeout: pollMs(30000)}).toBe(engine);
  await page.waitForFunction(() => !grok.shell.tv?.grid?.isRenderPending, null, {timeout: pollMs(30000)});
}, {description: 'Chem in this page runs the engine (its setting, read again within a second of a save), and the grid has repainted'});
