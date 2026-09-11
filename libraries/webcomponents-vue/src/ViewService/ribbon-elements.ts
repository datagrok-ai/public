import * as ui from 'datagrok-api/ui';
import {useDebounceFn} from '@vueuse/core';
import {DebounceFactory, DISABLED_STYLE_DEBOUNCE_TIME, RibbonMenuItem, RibbonPanelItem} from './ribbon-types';
import {dispatchClick, hoverText} from './ribbon-core';

export const defaultDebounce: DebounceFactory = (fn, wait) => useDebounceFn(fn, wait);

export interface ElementDeps {
  debounce: DebounceFactory;
  warn: (message: string) => void;
}

export interface PanelItemState {
  el: HTMLElement;
  iconKey: string;
  item: RibbonPanelItem;
  applyDisabledStyle: (disabled: boolean) => void;
}

export interface MenuItemState {
  el: HTMLElement;
  key: string;
  item: RibbonMenuItem;
}

export function createPanelItem(item: RibbonPanelItem, iconKey: string, deps: ElementDeps): PanelItemState {
  let el: HTMLElement;
  if (typeof item.icon === 'string')
    el = ui.iconFA(item.icon);
  else {
    const name = item.icon.path.split('/').pop()!.split('.')[0];
    el = ui.iconImage(name, item.icon.path);
    if (item.icon.width != null)
      el.style.width = `${item.icon.width}px`;
    if (item.icon.height != null)
      el.style.height = `${item.icon.height}px`;
  }
  const state: PanelItemState = {
    el,
    iconKey,
    item,
    applyDisabledStyle: deps.debounce((disabled: boolean) => {
      el.style.opacity = disabled ? '0.4' : '';
      el.style.filter = disabled ? 'grayscale(1)' : '';
    }, () => state.item.disabledStyleDebounce ?? DISABLED_STYLE_DEBOUNCE_TIME),
  };
  ui.tooltip.bind(el, () => hoverText(state.item));
  el.addEventListener('click', () => dispatchClick(state.item, deps.warn));
  return state;
}

export function updatePanelItem(state: PanelItemState, item: RibbonPanelItem): void {
  state.item = item;
  state.el.style.backgroundColor = item.active ? 'var(--grey-1)' : '';
  state.applyDisabledStyle(!!item.disabled && (item.disabledStyle ?? 'default') === 'default');
}

export function createMenuItem(item: RibbonMenuItem, key: string): MenuItemState {
  const el = document.createElement('div');
  if (item.icon) {
    const icon = ui.iconFA(item.icon);
    icon.style.width = '15px';
    icon.style.display = 'inline-block';
    icon.style.textAlign = 'center';
    el.append(icon, ' ');
  }
  el.append(item.text);
  const state: MenuItemState = {el, key, item};
  ui.tooltip.bind(el, () => hoverText(state.item));
  return state;
}
