import * as ui from 'datagrok-api/ui';

import {APP_NAME} from '../../common/view/const';
import {IsolatedAppUIBase} from '../../common/view/isolated-app-ui';
import {DataManager} from '../model/data-manager';
import {EventBus} from '../model/event-bus';
import {URLRouter} from '../model/router';
import {PatternAppLeftSection} from './components/left-section';
import {PatternAppRightSection} from './components/right-section';
import {PatternConfiguration} from '../model/types';


export class OligoPatternUI extends IsolatedAppUIBase {
  constructor() {
    super(APP_NAME.PATTERN);
  }

  protected getContent(): Promise<HTMLDivElement> {
    return getContent();
  }
}


async function getContent(): Promise<HTMLDivElement> {
  const dataManager = await DataManager.getInstance();
  const urlRouter = new URLRouter();

  const initialPatternConfig = await getInitialPatternConfig(dataManager, urlRouter);
  const eventBus = new EventBus(initialPatternConfig);
  urlRouter.updateURLOnPatternLoaded(eventBus);

  const leftSection = new PatternAppLeftSection(eventBus, dataManager).getLayout();
  const rightSection = new PatternAppRightSection(eventBus, dataManager).getLayout();

  const isResizeable = true;

  const layout = ui.splitH([
    leftSection,
    rightSection,
  ], {}, isResizeable);

  return layout;
}

async function getInitialPatternConfig(
  dataManager: DataManager,
  urlRouter: URLRouter
): Promise<PatternConfiguration> {
  const patternHash = urlRouter.getPatternHash();
  let initialPatternConfig = await dataManager.getPatternConfig(patternHash);
  if (!initialPatternConfig) {
    urlRouter.clearPatternURL();
    initialPatternConfig = dataManager.getDefaultPattern();
  }
  return initialPatternConfig;
}
