import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {RibbonHost} from './ribbon-types';
import {createRibbonService, RibbonService} from './ribbon-service';

export * from './ribbon-types';
export {createRibbonService} from './ribbon-service';
export type {RibbonService} from './ribbon-service';

export interface DgViewService extends RibbonService {
  readonly view: DG.ViewBase;
}

export const DG_VIEW_SERVICE_KEY: Vue.InjectionKey<DgViewService> = Symbol('DgViewService');

function dgRibbonHost(view: DG.ViewBase): RibbonHost {
  return {
    getPanels: () => view.getRibbonPanels(),
    setPanels: (panels) => view.setRibbonPanels(panels),
    get closing() {
      return view.closing;
    },
    rebuildGroup: (name, rank, els, onClick, opts) => {
      const group = view.ribbonMenu.group(name, rank);
      group.clear();
      group.items(els, onClick, opts);
    },
    removeGroup: (name) => {
      const menu = view.ribbonMenu;
      menu.group(name).clear();
      menu.remove(name);
    },
    warn: (message) => grok.shell.warning(message),
  };
}

export function createDgViewService(view: DG.ViewBase): DgViewService {
  const rawView = Vue.markRaw(view);
  return {view: rawView, ...createRibbonService(dgRibbonHost(rawView))};
}

export function provideDgViewService(app: Vue.App<any>, view: DG.ViewBase): DgViewService {
  const service = createDgViewService(view);
  app.provide(DG_VIEW_SERVICE_KEY, service);
  return service;
}

export function useViewService(): DgViewService {
  const service = Vue.inject(DG_VIEW_SERVICE_KEY);
  if (!service)
    throw new Error('DgViewService is not provided, call provideDgViewService(app, view) on the Vue app');
  return service;
}

export function useDgView(): DG.ViewBase {
  return useViewService().view;
}
