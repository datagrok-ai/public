import {OligoPatternUI} from '../../pattern/view/ui';
import {OligoStructureUI} from '../../structure/view/ui';
import {OligoTranslatorUI} from '../../translator/view/ui';
import {IsolatedAppUIBase} from './isolated-app-ui';
import {APP_NAME} from './const';

/** For plugins from external packages */
export class ExternalPluginUI extends IsolatedAppUIBase {
  constructor(viewName: string, private content: HTMLDivElement) {
    super(viewName);
  }

  protected getContent(): Promise<HTMLDivElement> {
    return Promise.resolve(this.content);
  }
}

export function getSpecifiedAppUI(appName: string): IsolatedAppUIBase {
  switch (appName) {
  case APP_NAME.TRANSLATOR:
    return new OligoTranslatorUI();
  case APP_NAME.PATTERN:
    return new OligoPatternUI();
  case APP_NAME.STRUCTURE:
    return new OligoStructureUI();
  default:
    throw new Error(`Unknown app name: ${appName}`);
  }
}
