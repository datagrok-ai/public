import * as Vue from 'vue';
import {DebounceFactory, RibbonHost, RibbonMenuItem, RibbonPanelItem}
  from '@datagrok-libraries/webcomponents-vue/src/ViewService/ribbon-types';

export interface HostCall {
  op: string;
  args: any[];
}

export interface FakeHost {
  host: RibbonHost;
  calls: HostCall[];
  setClosing(value: boolean): void;
}

export function fakeHost(initialPanels: HTMLElement[][] = []): FakeHost {
  const calls: HostCall[] = [];
  let closing = false;
  const host: RibbonHost = {
    get closing() {
      return closing;
    },
    getPanels: () => initialPanels,
    setPanels: (panels) => calls.push({op: 'setPanels', args: [panels]}),
    rebuildGroup: (name, rank, els, onClick, opts) =>
      calls.push({op: 'rebuildGroup', args: [name, rank, els, onClick, opts]}),
    removeGroup: (name) => calls.push({op: 'removeGroup', args: [name]}),
    warn: (message) => calls.push({op: 'warn', args: [message]}),
  };
  return {host, calls, setClosing: (value) => closing = value};
}

export const identityDebounce: DebounceFactory = (fn) => fn;

export const runInScope = (fn: () => void): Vue.EffectScope => {
  const scope = Vue.effectScope();
  scope.run(fn);
  return scope;
};

export const panelItem = (overrides: Partial<RibbonPanelItem> = {}): RibbonPanelItem =>
  ({icon: 'save', onClick: () => {}, ...overrides});

export const menuItem = (overrides: Partial<RibbonMenuItem> = {}): RibbonMenuItem =>
  ({text: 'item', onClick: () => {}, ...overrides});

export const countOp = (calls: HostCall[], op: string) => calls.filter((c) => c.op === op).length;

export const lastOp = (calls: HostCall[], op: string): HostCall => {
  const found = [...calls].reverse().find((c) => c.op === op);
  if (!found)
    throw new Error(`no recorded host call ${op}`);
  return found;
};
