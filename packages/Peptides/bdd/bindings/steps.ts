import {expect, type Page} from '@playwright/test';
import {Given, Then, When} from '@datagrok-libraries/bdd';
import {expectCustomEvent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

export const peptidesInitialized = Given('the Peptides package is initialized', async (page: Page) => {
  await page.evaluate(async () => { await grok.functions.call('Peptides:initPeptides'); });
}, {tier: 'api', description: 'waits for the sequence helper and monomer library initialization'});

export const sarReady = Then('the SAR analysis should be ready', async (page: Page) => {
  const event = await expectCustomEvent(page, 'peptides-sar-ready', 180000) as {table: string};
  expect(event.table).toBe(await page.evaluate(() => grok.shell.t.name));
}, {description: 'consumes peptides-sar-ready after launch or settings application, including MCL completion'});

export const sarSetting = Then('the SAR setting {string} should be {string}',
  async (page: Page, path: string, expected: string) => {
    const value = await page.evaluate((p) => {
      const settings = JSON.parse(grok.shell.t.getTag('settings'));
      return p.split('.').reduce((value, key) => value?.[key], settings);
    }, path);
    expect(String(value), `persisted SAR setting ${path}`).toBe(expected);
  }, {description: 'reads the persisted settings tag; pair with the resulting viewer or data change'});

export const openLanding = Given('user opens the Peptides landing view', async (page: Page) => {
  await page.evaluate(async () => {
    const view = await grok.functions.call('Peptides:Peptides');
    grok.shell.addView(view);
  });
}, {tier: 'api', description: 'opens the view returned by the Peptides entry function'});

export const openDemo = When('user opens the Peptide SAR demo dashboard', async (page: Page) => {
  await page.evaluate(async () => { await grok.functions.call('Peptides:macromoleculeSarFastaDemo'); });
}, {tier: 'api', description: 'runs the registered Peptide SAR dashboard entry point'});

export const sequenceSelection = Then('only rows with {string} at position {int} of {string} column should be selected',
  async (page: Page, monomer: string, position: number, column: string) => {
    const {matches, selected} = await page.evaluate(({monomer, position, column}) => {
      const df = grok.shell.t;
      const col = df.getCol(column);
      if (col.getTag('units') !== 'separator' || !col.getTag('separator'))
        throw new Error(`Column ${column} must use separator notation`);
      const separator = col.getTag('separator');
      const matches: number[] = [];
      const selected: number[] = [];
      for (let row = 0; row < df.rowCount; row++) {
        if (String(col.get(row)).split(separator)[position - 1] === monomer)
          matches.push(row);
        if (df.selection.get(row))
          selected.push(row);
      }
      return {matches, selected};
    }, {monomer, position, column});
    expect(matches.length, `${monomer} at position ${position} has matches`).toBeGreaterThan(0);
    expect(selected, 'the actual selected row indices').toEqual(matches);
  }, {description: 'compares the selection against every source separator sequence, before SAR adds position columns'});
