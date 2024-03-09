/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PatternAppLeftSection} from './left-section';
import {PatternAppRightSection} from './right-section';
import {PatternAppDataManager} from '../model/external-data-manager';
import {EventBus} from '../model/event-bus';
import {PatternDefaultsProvider} from '../model/defaults-provider';

export class PatternAppLayout {
  private leftSection: PatternAppLeftSection;
  private rightSection: PatternAppRightSection;

  constructor() {
    const defaultStateConfigurator = new PatternDefaultsProvider();
    const eventBus = new EventBus(defaultStateConfigurator);
    const dataManager = new PatternAppDataManager(eventBus);

    this.leftSection = new PatternAppLeftSection(eventBus, dataManager, defaultStateConfigurator);
    this.rightSection = new PatternAppRightSection(eventBus);
  }

  generateHTML(): HTMLDivElement {
    const leftSection = this.leftSection.getLayout();
    const rightSection = this.rightSection.getLayout();

    const isResizeable = true;

    const layout = ui.splitH([
      leftSection,
      rightSection,
    ], {}, isResizeable);

    return layout;
  }
}
