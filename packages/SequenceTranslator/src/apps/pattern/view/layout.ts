/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PatternAppLeftSection} from './left-section';
import {PatternAppRightSection} from './right-section';
import {PatternAppDataManager} from '../model/external-data-manager';
import {EventBus} from '../model/event-bus';
import {DefaultStateConfigurator} from '../model/default-state-configurator';

export class PatternAppLayout {
  private defaultStateConfigurator = new DefaultStateConfigurator();
  private eventBus = new EventBus(this.defaultStateConfigurator);
  private dataManager = new PatternAppDataManager(this.eventBus);
  private leftSection = new PatternAppLeftSection(this.eventBus, this.dataManager, this.defaultStateConfigurator);
  private rightSection = new PatternAppRightSection(this.eventBus, this.dataManager);

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
