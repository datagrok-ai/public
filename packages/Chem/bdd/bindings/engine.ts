/* Chem's Chemical Engine (RDKit or Crux, a package setting the stand's runs switch: misc/crux/engine-setting.sh), for
   the panes whose method is the engine's: Toxicity and Drug Likeness are OpenChemLib's in RDKit mode, RDKit's
   FilterCatalog alerts and QED in Crux mode (widgets/toxicity.ts, widgets/drug-likeness.ts). A journey runs its
   scenarios in one test, so a scenario cannot be skipped in one engine: its step checks what the engine shows. */
import {Page, test} from '@playwright/test';
import {Given, Then} from '@datagrok-libraries/bdd';
import {el, expect, locate, silent} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

/** The texts each engine's pane shows. */
const METHODS: {[engine: string]: {[pane: string]: string[]}} = {
  RDKit: {'Toxicity': ['Mutagenicity', 'Tumorigenicity', 'Irritating effects', 'Reproductive effects'],
    'Drug Likeness': ['Score:']},
  Crux: {'Toxicity': ['Brenk', 'NIH'], 'Drug Likeness': ['QED:', 'Aromatic rings']},
};

export const engineGate = Given('the Chemical Engine is {string}', async (page: Page, engine: string) => {
  silent(page);
  test.skip(await page.evaluate(() => grok.functions.call('Chem:chemicalEngine')) !== engine,
    `the Chemical Engine is not ${engine} — the rest of this test would run another engine's path`);
}, {tier: 'api', description: 'a gate: skips the rest of the test in the other engine; panes that are server-side Python in RDKit mode stay out of bdd'});

export const engineMethod = Then('{string} pane in {string} pane in context panel should show the Chemical Engine\'s method',
  async (page: Page, pane: string, group: string) => {
    const engine: string = await page.evaluate(() => grok.functions.call('Chem:chemicalEngine'));
    const loc = await locate(page, el(`"${pane}" pane in "${group}" pane in context panel`));
    for (const text of METHODS[engine]?.[pane] ?? [])
      await expect(loc).toContainText(text);
  }, {description: 'the pane shows its Chemical Engine\'s method: Toxicity OpenChemLib\'s four risks (RDKit) or the ' +
    'Brenk and NIH alerts (Crux); Drug Likeness OpenChemLib\'s score (RDKit) or QED and its properties (Crux)'});
