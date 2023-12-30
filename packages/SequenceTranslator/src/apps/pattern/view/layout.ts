/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PatternAppLeftSection} from './left-section';
import {PatternAppDataManager} from '../model/external-data-manager';
import {EventBus} from '../model/event-bus';
import {PatternConfigurationManager} from '../model/pattern-state-manager';

export class PatternLayoutController {
  private eventBus = new EventBus();
  private dataManager = PatternAppDataManager.getInstance(this.eventBus);
  private patternConfiguration = new PatternConfigurationManager(this.eventBus, this.dataManager);

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
    const leftSection = new PatternAppLeftSection(this.eventBus, this.dataManager, this.patternConfiguration);
    return leftSection.getLayout();
  }

  // private getRightSection(): HTMLDivElement {
  // }
}
