import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';


export namespace u2 {

  export function appHeader(iconPath: string, description: HTMLElement): HTMLElement {
    const icon = ui.iconImage('', iconPath);
    $(icon).addClass('ui-app-header-icon').css('margin-right', '20px');

    return ui.divH([
      icon,
      description
    ]);
  }

}