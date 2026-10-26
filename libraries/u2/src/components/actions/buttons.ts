/* Icon-capable buttons. They live here rather than in core/elements.ts because icons are a
   component (L2) and core may not import upwards — both factories compose the core `button()`,
   so the plain and the icon-bearing button stay one skin. */
import {Signal} from '../../core/signals.js';
import {Scope, debug} from '../../core/scope.js';
import {Tooltip} from '../../core/tooltip.js';
import {button, div, span, Text} from '../../core/elements.js';
import {icon, IconVariant} from '../display/icon.js';
import {Menu} from '../navigation/menu.js';

export interface IconButtonOptions {
  variant?: IconVariant;
  tooltip?: string;
  /** Pressed state: the click writes it and an effect mirrors it into `aria-pressed`, so a
   * toggling button needs an ambient Scope. */
  toggle?: Signal<boolean>;
  /** A menu built on every press and shown under the button; pressed again, it closes. */
  menu?: () => Menu;
}

const NO_OWNER = 'u2: signal binding needs an owner — ' +
  'wrap the code in Control.build(...) or component.runInScope(...)';

const warned = new Set<string>();

/** Builds under a scope of the factory's own: a child of the ambient one, so it goes with whoever
 * built it, releases its listeners and popups together and is named in leak reports; none outside
 * any scope (the element's own listeners then live as long as the element). */
export function ownScope<T>(label: string, build: (scope: Scope | undefined) => T): T {
  const ambient = Scope.ambient;
  const scope = ambient === undefined ? undefined : new Scope(ambient, {label, lazy: true});
  return Scope.runWith(scope, () => build(scope));
}

/** Hover text in the icon.ts convention: the tooltip service when there is a scope to own it,
 * the native `title` otherwise. The tooltip doubles as the accessible name — the u2 tooltip
 * leaves no `title` attribute, so without `aria-label` an icon-only button would be nameless. */
function tip(el: HTMLElement, text: string): void {
  el.setAttribute('aria-label', text);
  const scope = Scope.ambient;
  if (scope)
    Tooltip.bind(el, text, scope);
  else
    el.title = text;
}

/** Square icon-only button — flat until hover, sized by its context (the form metric standalone,
 * the inline metric inside a toolbar or a flat button group). */
export function iconButton(name: string, onClick: (e: MouseEvent) => void,
  options?: IconButtonOptions): HTMLButtonElement {
  const toggle = options?.toggle;
  let menu: Menu | undefined;
  const el = button(icon(name, {variant: options?.variant}), (e) => {
    if (toggle)
      toggle.value = !toggle.peek();
    if (options?.menu !== undefined) {
      const open = menu?.isOpen.peek() === true;
      menu?.close();
      menu = open ? undefined : options.menu();
      menu?.show({anchor: el});
    }
    onClick(e);
  });
  el.classList.add('u2-icon-btn');
  el.dataset.u2 = 'icon-button';
  if (options?.tooltip)
    tip(el, options.tooltip);
  else if (debug.enabled && !warned.has(name)) {
    warned.add(name);
    console.warn(`u2: iconButton("${name}") has no tooltip — an icon-only button needs one for its accessible name`);
  }
  if (!toggle)
    return el;
  const scope = Scope.ambient;
  if (!scope)
    throw new Error(NO_OWNER);
  scope.effect(() => {
    const pressed = toggle.value;
    el.classList.toggle('u2-icon-btn-pressed', pressed);
    el.setAttribute('aria-pressed', String(pressed));
  });
  return el;
}

/** Text button with a leading (or trailing) icon. The label is its own span so a signal-bound
 * caption rewrites the label alone and never the icon. */
export function buttonWithIcon(text: Text, iconName: string, onClick: (e: MouseEvent) => void,
  options?: {primary?: boolean, iconRight?: boolean, tooltip?: string}): HTMLButtonElement {
  const el = button(icon(iconName), onClick, {primary: options?.primary});
  const label = span(text, 'u2-btn-label');
  if (options?.iconRight)
    el.prepend(label);
  else
    el.append(label);
  el.dataset.u2 = 'button';
  if (options?.tooltip)
    tip(el, options.tooltip);
  return el;
}

export interface DropDownButtonOptions {
  icon?: string;
  primary?: boolean;
  tooltip?: string;
  /** Split mode: the main part fires this, and only the arrow segment opens the menu. */
  split?: (e: MouseEvent) => void;
}

/** Button that opens a menu, or — with `split` — a main action plus an attached arrow that opens
 * it. `buildMenu` runs on every open, so checked items never show stale state. Built under a scope
 * of its own (see {@link ownScope}): its listeners and an open menu go together.
 * `aria-expanded` follows the menu through a per-open Scope holding one effect on `isOpen`, which
 * that same effect disposes on close, so nothing outlives the menu. */
export function dropDownButton(text: Text, buildMenu: (menu: Menu) => void,
  options?: DropDownButtonOptions): HTMLElement {
  return ownScope('dropdown-button', (scope) => {
    const split = options?.split;
    const primary = options?.primary;
    let root: HTMLElement;
    let trigger: HTMLElement;
    let current: Menu | undefined;

    const close = (): boolean => {
      if (!current)
        return false;
      current.close();
      return true;
    };

    const toggle = (): void => {
      if (close())
        return;
      const menu = Menu.popup();
      buildMenu(menu);
      current = menu;
      menu.show({anchor: root});
      const watch = new Scope(undefined, {label: 'dropdown-menu'});
      watch.effect(() => {
        const open = menu.isOpen.value;
        trigger.setAttribute('aria-expanded', String(open));
        if (open)
          return;
        if (current === menu)
          current = undefined;
        watch.dispose();
      });
    };

    const onMain = split ? (e: MouseEvent) => {
      close();
      split(e);
    } : toggle;
    const main = options?.icon ?
      buttonWithIcon(text, options.icon, onMain, {primary}) :
      button(span(text, 'u2-btn-label'), onMain, {primary});
    if (options?.tooltip)
      tip(main, options.tooltip);

    if (split) {
      const arrow = button(icon('chevron-down', {cls: 'u2-dropdown-chevron'}), toggle, {primary});
      arrow.classList.add('u2-dropdown-arrow');
      arrow.setAttribute('aria-label', 'More');
      root = div([main, arrow], 'u2-dropdown-split');
      trigger = arrow;
    } else {
      main.append(icon('chevron-down', {cls: 'u2-dropdown-chevron'}));
      main.classList.add('u2-dropdown-btn');
      root = main;
      trigger = main;
    }
    root.dataset.u2 = 'dropdown-button';
    trigger.setAttribute('aria-haspopup', 'menu');
    trigger.setAttribute('aria-expanded', 'false');

    const onKeyDown = (e: Event) => {
      if ((e as KeyboardEvent).key !== 'ArrowDown' || current)
        return;
      e.preventDefault();
      toggle();
    };
    if (scope) {
      scope.listen(root, 'keydown', onKeyDown);
      scope.own(close);
    } else
      root.addEventListener('keydown', onKeyDown);
    return root;
  });
}
