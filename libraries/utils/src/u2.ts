import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import $ from 'cash-dom';
import './u2.css';


export namespace u2 {

  /** */
  export namespace panels {

    /** Vertically positioned items */
    export function vert(items: any[]): HTMLElement {
      return ui.divV(items.map((item) => ui.render(item)), { classes: 'u2-panel'});
    }

    /** Horizontally positioned items */
    export function horz(items: any[]): HTMLElement {
      return ui.divH(items.map((item) => ui.render(item)), { classes: 'u2-panel'});
    }

  }

  export interface IAppInfo {
    iconPath: string;
    description: string;
    learnMoreUrl?: string;
    appTitle?: string;
    appSubTitle?: string;
    bottomLine?: boolean;
    iconSize?: number;
  }

  export function appHeader(header: IAppInfo): HTMLElement {
    const icon = ui.iconImage('', header.iconPath);
    $(icon).addClass('ui-app-header-icon').css('margin-right', '20px');
    if (header.iconSize) {
      icon.style.width = `${header.iconSize}px`;
      icon.style.height = `${header.iconSize}px`;
    }
    const appHeaderDiv = panels.horz([]);
    if (header.appTitle) {
      appHeaderDiv.classList.add('u2-app-header-with-name');
      const nameDiv = ui.divV([ui.divText(header.appTitle, 'u2-app-header-app-name-div')], 'u2-app-header-app-name-and-slogan-div');
      if (header.appSubTitle)
        nameDiv.append(ui.divText(header.appSubTitle));
      const iconDiv = ui.div(icon);
      appHeaderDiv.append(ui.divH([iconDiv, nameDiv]));
    }
    else
      appHeaderDiv.append(icon);

    const description = panels.vert([
      ui.markdown(header.description),
      header.learnMoreUrl ? ui.link('Learn more', header.learnMoreUrl, undefined, {style: {marginLeft: '15px'}}) : null
    ]);
    description.classList.add('u2-app-header-description');
    appHeaderDiv.append(description);
    if (header.bottomLine)
      appHeaderDiv.classList.add('u2-app-header-bottom-line');
    return appHeaderDiv;
  }


  export namespace tools {

    /** Executes {@link func} while showing the "running" indicator on {@link root}.
     * Handles and logs exceptions. */
    export async function runAsync<T>(root: HTMLElement, func: () => Promise<T>) {
      ui.setUpdateIndicator(root, true);
      try {
        return await func();
      }
      catch (e) {
        grok.log.error(e as string);
      }
      finally {
        ui.setUpdateIndicator(root, false);
      }
    }
  }
}