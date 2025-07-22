import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export class HelpObjectHandler extends DG.ObjectHandler {
  isApplicable(x: any): boolean {
    return x instanceof HelpObject;
  }

  get type() {
    return 'PSHelpObject';
  }

  renderMarkup(x: HelpObject, context?: any): HTMLElement {
    const icon = ui.icons.info(() => {});
    const title = ui.label(x.title);
    const container = ui.divH([icon, title], {style: {alignItems: 'center'}});
    container.addEventListener('click', (e) => {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
      grok.shell.windows.showHelp = true;
      grok.shell.windows.help.showHelp(x.helpURL);
    });

    container.addEventListener('dblclick', (e) => {
      e.preventDefault();
      e.stopPropagation();
      e.stopImmediatePropagation();
      window.open('https://datagrok.ai' + x.helpURL, '_blank');
    });
    return container;
  }

  renderTooltip(x: HelpObject, context?: any): HTMLElement {
    return ui.divV([
      ui.h3(x.helpURL.split('/').slice(2, -1).join(' | ') + ' | ' + x.title, {style: {marginTop: '0px'}}),
      ui.divText('Click to open help panel'),
      ui.divText('Double-click to open in a new tab')]);
  }
}

export class HelpObject {
  public readonly title: string;
  public readonly helpURL: string;
  public readonly keywords?: string[];

  constructor(title: string, helpURL: string, keywords?: string[]) {
    this.title = title;
    this.helpURL = helpURL;
    this.keywords = keywords;
  }

  static fromHelpInfo(helpInfo: {title: string, helpURL: string, keywords?: string[]}): HelpObject {
    return new HelpObject(helpInfo.title, helpInfo.helpURL, helpInfo.keywords);
  }
}
