import {RibbonItemBase, RibbonMenuItem, RibbonPanelItem} from './ribbon-types';

export const resolveReason = (item: RibbonItemBase): string | null =>
  (typeof item.disabledReason === 'function' ? item.disabledReason() : item.disabledReason) ?? null;

export type ClickResolution = {kind: 'run'} | {kind: 'warn', reason: string} | {kind: 'ignore'};

export const resolveClick = (item: RibbonItemBase): ClickResolution => {
  if (!item.disabled)
    return {kind: 'run'};
  const mode = item.disabledReasonMode ?? 'tooltip';
  if (mode === 'popup' || mode === 'both') {
    const reason = resolveReason(item);
    if (reason)
      return {kind: 'warn', reason};
  }
  return {kind: 'ignore'};
};

export const dispatchClick = (item: RibbonItemBase, warn: (reason: string) => void): void => {
  const resolution = resolveClick(item);
  if (resolution.kind === 'run')
    item.onClick();
  else if (resolution.kind === 'warn')
    warn(resolution.reason);
};

export const tooltipReason = (item: RibbonItemBase): string | null => {
  const mode = item.disabledReasonMode ?? 'tooltip';
  return item.disabled && (mode === 'tooltip' || mode === 'both') ? resolveReason(item) : null;
};

export const hoverText = (item: RibbonItemBase): string | null =>
  tooltipReason(item) ?? item.tooltip ?? null;

export const bySortKey = (a: {priority?: number, seq: number}, b: {priority?: number, seq: number}) =>
  ((a.priority ?? Infinity) - (b.priority ?? Infinity)) || (a.seq - b.seq);

export const panelIconKey = (icon: RibbonPanelItem['icon']) =>
  typeof icon === 'string' ? `fa:${icon}` : `img:${icon.path}`;

export const menuItemKey = (item: RibbonMenuItem) => `${item.text}|${item.icon ?? ''}`;

export interface MenuRegSnapshot {
  seq: number;
  priority?: number;
  name: string;
  items: RibbonMenuItem[];
}

// Merges same-name contributions in registration order; group rank comes from
// the first-seen contributor; empty groups are dropped.
export function composeMenuGroups(regs: MenuRegSnapshot[]): {name: string, items: RibbonMenuItem[]}[] {
  const byName = new Map<string, {priority?: number, seq: number, items: RibbonMenuItem[]}>();
  for (const reg of [...regs].sort((a, b) => a.seq - b.seq)) {
    const group = byName.get(reg.name);
    if (group)
      group.items.push(...reg.items);
    else
      byName.set(reg.name, {priority: reg.priority, seq: reg.seq, items: [...reg.items]});
  }
  return [...byName.entries()]
    .filter(([, group]) => group.items.length > 0)
    .sort(([, a], [, b]) => bySortKey(a, b))
    .map(([name, group]) => ({name, items: group.items}));
}

// Keyed reconcile: reuses prev states by key (duplicates pooled in order),
// creates the rest; `changed` reflects any creation, removal, or reorder.
export function reconcileByKey<S, I>(
  prev: S[],
  next: I[],
  keyOf: (item: I) => string,
  keyOfState: (state: S) => string,
  create: (item: I, key: string) => S,
  update: (state: S, item: I) => void,
): {states: S[], changed: boolean} {
  const pools = new Map<string, S[]>();
  for (const state of prev) {
    const key = keyOfState(state);
    const pool = pools.get(key);
    if (pool)
      pool.push(state);
    else
      pools.set(key, [state]);
  }
  const states = next.map((item) => {
    const key = keyOf(item);
    const state = pools.get(key)?.shift() ?? create(item, key);
    update(state, item);
    return state;
  });
  const changed = states.length !== prev.length || states.some((state, i) => state !== prev[i]);
  return {states, changed};
}

export const panelsEqual = (a: HTMLElement[][], b: HTMLElement[][]) =>
  a.length === b.length &&
  a.every((panel, i) => panel.length === b[i].length && panel.every((el, j) => el === b[i][j]));
