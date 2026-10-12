/* The monomer manager's editor form (Bio > Manage > Monomers): its R-groups tab holds a grid of the R groups the
   drawn monomer has (an ItemsGrid: a row of inputs per R group, its cap's SMILES, alternate ID, cap name and label),
   filled from the R labels of the sketcher's SMILES. */
import type {Page} from '@playwright/test';
import {element, Then} from '@datagrok-libraries/bdd';
import {expect} from '@datagrok-libraries/bdd/runtime';

element('monomer sketcher', {selector: '.monomer-manager-sketcher',
  description: 'the editor\'s sketcher host (DG.chem.Sketcher): its molecule field, its options menu (≡) and the sketcher in it; made once, with the sketcher of that moment, and kept for the page\'s life'});

export const rgroupLabels = Then('the R-groups grid of the monomer form should list {string}', async (page: Page, list: string) => {
  const want = list.split(/\s*,\s*/).filter(Boolean);
  await expect.poll(() => page.evaluate(() => {
    const grid = document.querySelector('.monomer-manager-form-tab-control .ui-items-grid');
    if (!grid)
      return ['no R-groups grid'];
    // a row's label is the input that holds R and a number (its alternate ID is R1-H, its cap SMILES [*:1][H])
    return Array.from(grid.querySelectorAll('input')).map((i) => (i as HTMLInputElement).value)
      .filter((v) => /^R\d+$/.test(v));
  }), {message: 'the labels of the R-groups grid'}).toEqual(want);
}, {description: 'the labels of the grid\'s rows, in order (the inputs that hold R and a number)'});
