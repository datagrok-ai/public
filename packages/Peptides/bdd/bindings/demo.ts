import {expect, type Page} from '@playwright/test';
import {Then} from '@datagrok-libraries/bdd';

declare const DG: any;

export const demoRegistration = Then('the Peptide SAR demo should be registered as a Bioinformatics dashboard',
  async (page: Page) => {
    const options = await page.evaluate(() => {
      const func = DG.Func.find({package: 'Peptides', name: 'macromoleculeSarFastaDemo'})[0];
      return {path: func?.options.demoPath, dashboard: func?.options.isDemoDashboard};
    });
    expect(options).toEqual({path: 'Bioinformatics | Peptide SAR', dashboard: 'true'});
  }, {description: 'the registered function has the gallery path Bioinformatics | Peptide SAR and isDemoDashboard=true'});
