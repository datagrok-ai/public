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
  static async generateHTML(): Promise<HTMLDivElement> {
    const defaultsProvider = new PatternDefaultsProvider();
    const eventBus = new EventBus(defaultsProvider);

    const dataManager = await PatternAppDataManager.getInstance(eventBus);

    const leftSection = new PatternAppLeftSection(eventBus, dataManager, defaultsProvider).getLayout();
    const rightSection = new PatternAppRightSection(eventBus, dataManager).getLayout();

    const isResizeable = true;

    const layout = ui.splitH([
      leftSection,
      rightSection,
    ], {}, isResizeable);

    return layout;
  }
}
