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
  }

  export function appHeader(header: IAppInfo): HTMLElement {
    const icon = ui.iconImage('', header.iconPath);
    $(icon).addClass('ui-app-header-icon').css('margin-right', '20px');

    return panels.horz([
      icon,
      panels.vert([
        ui.markdown(header.description),
        header.learnMoreUrl ? ui.link('Learn more', header.learnMoreUrl) : null
      ])
    ]);
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