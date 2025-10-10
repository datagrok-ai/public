import * as ui from 'datagrok-api/ui';
import $ from 'cash-dom';

export function prop<T>(
  valueFunc: () => T,
  host: HTMLElement | null,
  tooltip: string,
  onClick?: () => Promise<void> | void
): { value: T, addColumnIcon: HTMLElement | null } {
  let addColumnIcon: HTMLElement | null = null;

  if (host && onClick) {
    addColumnIcon = ui.iconFA('plus', onClick, tooltip);

    ui.tools.setHoverVisibility(host, [addColumnIcon]);
    $(addColumnIcon)
      .css('color', '#2083d5')
      .css('position', 'absolute')
      .css('top', '2px')
      .css('left', '-12px')
      .css('margin-right', '5px');
  }

  return { value: valueFunc(), addColumnIcon };
}
