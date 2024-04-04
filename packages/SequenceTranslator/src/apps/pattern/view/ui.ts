import * as ui from 'datagrok-api/ui';

import {APP_NAME} from '../../common/view/const';
import {IsolatedAppUIBase} from '../../common/view/isolated-app-ui';
import {PatternConfigManager} from '../model/config-manager';
import {DataInitializer} from '../model/data-initializer';
import {PatternDefaultsProvider} from '../model/defaults-provider';
import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {URLRouter} from '../model/router';
import {PatternAppLeftSection} from './components/left-section';
import {PatternAppRightSection} from './components/right-section';

class PatternApp {
  static async getContent(): Promise<HTMLDivElement> {
    const defaultsProvider = new PatternDefaultsProvider();
    const dataInitializer = await DataInitializer.getInstance(defaultsProvider);
    const patternConfigManager = new PatternConfigManager(dataInitializer);

    const searchParams = new URLSearchParams(window.location.search);

    const eventBus = new EventBus(defaultsProvider);
    const dataManager = new PatternAppDataManager(eventBus, patternConfigManager);

    const urlRouter = new URLRouter(eventBus, dataManager);
    await urlRouter.navigate();

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

export class OligoPatternUI extends IsolatedAppUIBase {
  constructor() {
    super(APP_NAME.PATTERN);
  }

  protected getContent(): Promise<HTMLDivElement> {
    return PatternApp.getContent();
  }
}

