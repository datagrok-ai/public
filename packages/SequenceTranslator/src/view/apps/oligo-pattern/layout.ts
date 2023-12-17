/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {LeftSection} from './left-section';
import {DataManager} from './utils';

export class PatternLayoutHandler {
  private dataManager = new DataManager();

  get htmlDivElement(): HTMLDivElement {
    const leftSection = this.getLeftSection();
    // const rightSection = this.getRightSection();
    const layout = ui.splitH([
      leftSection,
      // rightSection
    ], {}, true);

    return layout;
  }

  private getLeftSection(): HTMLDivElement {
    return new LeftSection(this.dataManager).getLayout();
  }

  // private getRightSection(): HTMLDivElement {
  // }
}
