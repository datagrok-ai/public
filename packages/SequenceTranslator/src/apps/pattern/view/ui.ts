import * as ui from 'datagrok-api/ui';

import {PatternDefaultsProvider} from '../model/defaults-provider';
import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternAppLeftSection} from './components/left-section';
import {PatternAppRightSection} from './components/right-section';

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
