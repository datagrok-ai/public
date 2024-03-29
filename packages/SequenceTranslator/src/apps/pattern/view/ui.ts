import * as ui from 'datagrok-api/ui';

import {APP_NAME} from '../../common/view/const';
import {IsolatedAppUIBase} from '../../common/view/isolated-app-ui';
import {Router} from '../model/router';
import {PatternDefaultsProvider} from '../model/defaults-provider';
import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternAppLeftSection} from './components/left-section';
import {PatternAppRightSection} from './components/right-section';

class PatternApp {
  static async getContent(): Promise<HTMLDivElement> {
    const defaultsProvider = new PatternDefaultsProvider();
    const eventBus = new EventBus(defaultsProvider);
    const dataManager = await PatternAppDataManager.getInstance(eventBus);

    const urlRouter = new Router(eventBus, dataManager);
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

