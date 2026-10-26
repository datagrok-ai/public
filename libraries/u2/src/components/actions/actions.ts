/* Item actions — one declaration feeding every surface (docs/recipes/list-item-rendering.md):
   the hover-revealed icon block on a row shows the subset with icons, the right-click menu always
   shows the full list. */
import {iconButton, ownScope} from './buttons.js';
import {Menu} from '../navigation/menu.js';
import {Access} from '../../core/access.js';
import type {Capability} from '../../core/access.js';

export interface Action {
  name: string;
  /** Platform icon name (see `icon()`); an action without one appears only in menus. */
  icon?: string;
  /** An action that does not apply to the current object stays visible and greyed out. */
  enabled?: boolean;
  /** The capability the action needs (`'delete'`, a custom permission name): denied by the
   * access in force, the action is not rendered at all — permission ⇒ hidden, state ⇒ disabled. */
  requires?: Capability;
  /** The key shown beside it in menus ("Del"). */
  shortcut?: string;
  /** Grey text after the name in a menu ("3 issues"). */
  detail?: string;
  /** Throws something away: after a separator in menus, red in the row's icons. */
  destructive?: boolean;
  /** Starts a new section of a menu: a separator above it. */
  separator?: boolean;
  run: () => void;
}

export interface ActionsOptions {
  /** Default {@link Access.full}, so nothing outside EMS changes. */
  access?: Access;
  /** The object the actions apply to; refines `access` per row where the row carries its own. */
  row?: unknown;
}

/** Permission ⇒ hidden: the actions whose `requires` the access in force — refined by the row where
 * one is given — grants. What {@link rowActions} and {@link actionsMenu} render; a caller feeding a
 * list's `contextActions` applies it itself. */
export function allowedActions(actions: Action[], options?: ActionsOptions): Action[] {
  const table = options?.access ?? Access.full;
  const access = options?.row === undefined ? table : table.row(options.row);
  const out: Action[] = [];
  // a denied section head passes its separator to the next allowed action
  let separator = false;
  for (const action of actions) {
    if (action.requires !== undefined && !access.can(action.requires)) {
      separator ||= action.separator === true;
      continue;
    }
    out.push(separator && !action.separator ? {...action, separator: true} : action);
    separator = false;
  }
  return out;
}

function addItems(menu: Menu, actions: Action[]): void {
  for (const [i, action] of actions.entries()) {
    const previous = actions[i - 1];
    if (previous !== undefined && (action.separator || action.destructive && !previous.destructive))
      menu.separator();
    menu.item(action.name, action.run,
      {icon: action.icon, enabled: action.enabled, shortcut: action.shortcut, detail: action.detail});
  }
}

/** The hover-revealed action block, Slack/Outlook style: an icon button per icon-bearing action,
 * invisible until the enclosing `.u2-list-row` — or any element carrying `.u2-row` — is hovered,
 * or a button in it is tabbed to (css/list.css: opacity, not display, so the buttons keep their
 * space and stay focusable). Shortcuts only: the same list belongs on the right-click menu too. */
export function rowActions(actions: Action[], options?: ActionsOptions): HTMLElement {
  return ownScope('row-actions', () => {
    const el = document.createElement('div');
    el.className = 'u2-row-actions';
    el.dataset.u2 = 'row-actions';
    el.setAttribute('role', 'toolbar');
    for (const action of allowedActions(actions, options)) {
      if (!action.icon)
        continue;
      const button = iconButton(action.icon, () => action.run(), {tooltip: action.name});
      button.dataset.u2Item = action.name;
      if (action.destructive)
        button.classList.add('u2-row-action-destructive');
      if (action.enabled === false)
        button.disabled = true;
      el.append(button);
    }
    return el;
  });
}

/** A submenu under the flat list: what a wider target than the clicked object offers — the rows of
 * a multi-selection, run as one. */
export interface ActionGroup {
  label: string;
  actions: Action[];
}

export interface ActionsMenuOptions extends ActionsOptions {
  groups?: ActionGroup[];
}

/** A popup menu over the same action list — the right-click superset of {@link rowActions}. */
export function actionsMenu(actions: Action[], options?: ActionsMenuOptions): Menu {
  const menu = new Menu();
  const flat = allowedActions(actions, options);
  addItems(menu, flat);
  let anything = flat.length > 0;
  for (const group of options?.groups ?? []) {
    const allowed = allowedActions(group.actions, options);
    if (allowed.length === 0)
      continue;
    if (anything)
      menu.separator();
    anything = true;
    menu.group(group.label);
    addItems(menu, allowed);
    menu.endGroup();
  }
  return menu;
}
