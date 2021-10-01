import * as DG from 'datagrok-api/dg';
import { filter } from 'rxjs/operators';
import { Tutorial } from '../../../tutorial';


export class FiltersTutorial extends Tutorial {
  get name() {
    return 'Filters';
  }
  get description() {
    return 'A set of controls for quick filtering, selection, ' +
      'and visual assessment of column values';
  }

  get steps() { return 3; }
  
  protected async _run() {
    this.header.textContent = 'Filters';
    this.describe('Dynamic filtering is an important concept in exploratory data analysis, and ' +
      'our platform makes it as powerful and easy to use as possible. Let\'s start with opening ' +
      'a couple of regular viewers, so that effects of filtering would be immediately visible.');

    await this.openPlot('scatter plot', (x) => x.type === DG.VIEWER.SCATTER_PLOT);
    await this.openPlot('histogram', (x) => x.type === DG.VIEWER.HISTOGRAM);

    const filterDescription = `Now, let's add filters.\nWhile appearing very simple, this is a very ` +
      'powerful tool. It combines multiple indicators along with a few controls that let you easily ' +
      `modify current filters and selection. Let's start with filtering by a specific category.`;
    const filters = await this.openPlot('filters', (x) => x.type === DG.VIEWER.FILTERS, filterDescription);

    const filters = await this.openPlot('filters', `Now, let's add filters.\nWhile appearing very simple, this is a very powerful tool. It combines multiple indicators along with a few controls that let you easily modify current filters and selection. Let's start with filtering by a specific category.`, (x) => x.type === DG.VIEWER.FILTERS);

  }
}
