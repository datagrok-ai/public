import * as Vue from 'vue';
import {DebounceFactory, RibbonHost, RibbonMenuItem, RibbonPanelItem} from './ribbon-types';
import {bySortKey, composeMenuGroups, dispatchClick, menuItemKey, panelIconKey, panelsEqual,
  reconcileByKey, tooltipReason} from './ribbon-core';
import {createMenuItem, createPanelItem, defaultDebounce, MenuItemState,
  PanelItemState, updatePanelItem} from './ribbon-elements';

export interface RibbonService {
  registerPanel(items: () => RibbonPanelItem[], priority?: () => number | undefined): void;
  registerMenuGroup(name: () => string, items: () => RibbonMenuItem[], priority?: () => number | undefined): void;
  dispose(): void;
}

interface Registration<T> {
  seq: number;
  items: () => T[];
  priority?: () => number | undefined;
}

interface MenuRegistration extends Registration<RibbonMenuItem> {
  name: () => string;
}

export function createRibbonService(host: RibbonHost, opts?: {debounce?: DebounceFactory}): RibbonService {
  const deps = {debounce: opts?.debounce ?? defaultDebounce, warn: (message: string) => host.warn(message)};
  const scope = Vue.effectScope(true);
  const panelRegs = Vue.shallowReactive(new Map<symbol, Registration<RibbonPanelItem>>());
  const menuRegs = Vue.shallowReactive(new Map<symbol, MenuRegistration>());
  const externalPanels = host.getPanels();
  let seqCounter = 0;

  let panelStates = new Map<symbol, PanelItemState[]>();
  let lastPanels: HTMLElement[][] = [];

  let menuStates = new Map<string, MenuItemState[]>();
  const elToMenuState = new WeakMap<HTMLElement, MenuItemState>();
  let lastGroupNames: string[] = [];

  const menuDispatch = (el: HTMLElement) => {
    const state = elToMenuState.get(el);
    if (state)
      dispatchClick(state.item, deps.warn);
  };

  const menuIsValid = (el: HTMLElement): string | null => {
    const item = elToMenuState.get(el)?.item;
    if (!item?.disabled || (item.disabledStyle ?? 'default') === 'none')
      return null;
    return tooltipReason(item) ?? '';
  };

  const menuIsChecked = (el: HTMLElement) => !!elToMenuState.get(el)?.item.check;

  scope.run(() => {
    Vue.watch(
      () => [...panelRegs.entries()]
        .map(([token, reg]) => ({token, seq: reg.seq, priority: reg.priority?.(), items: reg.items()}))
        .sort(bySortKey),
      (next) => {
        const nextStates = new Map<symbol, PanelItemState[]>();
        const composed: HTMLElement[][] = [];
        for (const reg of next) {
          const {states} = reconcileByKey(
            panelStates.get(reg.token) ?? [], reg.items,
            (item) => panelIconKey(item.icon), (state) => state.iconKey,
            (item, key) => createPanelItem(item, key, deps), updatePanelItem);
          nextStates.set(reg.token, states);
          if (states.length)
            composed.push(states.map((s) => s.el));
        }
        panelStates = nextStates;
        if (host.closing)
          return;
        if (!panelsEqual(composed, lastPanels))
          host.setPanels([...externalPanels, ...composed]);
        lastPanels = composed;
      }, {flush: 'post'});

    Vue.watch(
      () => composeMenuGroups([...menuRegs.values()]
        .map((reg) => ({seq: reg.seq, priority: reg.priority?.(), name: reg.name(), items: reg.items()}))),
      (next) => {
        const nextStates = new Map<string, MenuItemState[]>();
        const changedGroups: string[] = [];
        for (const group of next) {
          const {states, changed} = reconcileByKey(
            menuStates.get(group.name) ?? [], group.items,
            menuItemKey, (state) => state.key,
            (item, key) => {
              const state = createMenuItem(item, key);
              elToMenuState.set(state.el, state);
              return state;
            },
            (state, item) => {
              state.item = item;
            });
          nextStates.set(group.name, states);
          if (changed)
            changedGroups.push(group.name);
        }
        menuStates = nextStates;
        const nextNames = next.map((g) => g.name);
        if (host.closing) {
          lastGroupNames = nextNames;
          return;
        }
        const orderChanged = nextNames.length !== lastGroupNames.length ||
          nextNames.some((name, i) => name !== lastGroupNames[i]);
        const rebuild = (name: string, rank: number) => {
          const states = nextStates.get(name)!;
          const anyCheck = states.some((s) => s.item.check !== undefined);
          host.rebuildGroup(name, rank, states.map((s) => s.el), menuDispatch,
            {isValid: menuIsValid, ...(anyCheck ? {isChecked: menuIsChecked} : {})});
        };
        if (orderChanged) {
          for (const name of lastGroupNames)
            host.removeGroup(name);
          next.forEach((group, rank) => rebuild(group.name, rank));
        } else {
          for (const name of changedGroups)
            rebuild(name, nextNames.indexOf(name));
        }
        lastGroupNames = nextNames;
      }, {flush: 'post'});
  });

  const register = <T>(registry: Map<symbol, T>, entry: T) => {
    const token = Symbol();
    registry.set(token, entry);
    if (Vue.getCurrentScope())
      Vue.onScopeDispose(() => registry.delete(token));
  };

  return {
    registerPanel: (items, priority) =>
      register(panelRegs, {seq: seqCounter++, items, priority}),
    registerMenuGroup: (name, items, priority) =>
      register(menuRegs, {seq: seqCounter++, name, items, priority}),
    dispose: () => {
      scope.stop();
      if (host.closing)
        return;
      host.setPanels(externalPanels);
      for (const name of lastGroupNames)
        host.removeGroup(name);
    },
  };
}
