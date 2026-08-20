/* Item actions — one declaration feeding every surface (docs/recipes/list-item-rendering.md):
   the hover-revealed icon block on a row shows the subset with icons, the right-click menu always
   shows the full list. */
import {iconButton} from './buttons.js';
import {Menu} from './menu.js';

export interface Action {
  name: string;
  /** Platform icon name (see `icon()`); an action without one appears only in menus. */
  icon?: string;
  /** An action that does not apply to the current object stays visible and greyed out. */
  enabled?: boolean;
  run: () => void;
}

/** The hover-revealed action block, Slack/Outlook style: an icon button per icon-bearing action,
 * invisible until the enclosing `.u2-list-row` — or any element carrying `.u2-row` — is hovered,
 * or a button in it is tabbed to (css/list.css: opacity, not display, so the buttons keep their
 * space and stay focusable). Shortcuts only: the same list belongs on the right-click menu too. */
export function rowActions(actions: Action[]): HTMLElement {
  const el = document.createElement('div');
  el.className = 'u2-row-actions';
  el.dataset.u2 = 'row-actions';
  for (const action of actions) {
    if (!action.icon)
      continue;
    const button = iconButton(action.icon, action.run, {tooltip: action.name});
    if (action.enabled === false)
      button.disabled = true;
    el.append(button);
  }
  return el;
}

/** A popup menu over the same action list — the right-click superset of {@link rowActions}. */
export function actionsMenu(actions: Action[]): Menu {
  const menu = new Menu();
  for (const action of actions)
    menu.item(action.name, action.run, {icon: action.icon, enabled: action.enabled});
  return menu;
}
