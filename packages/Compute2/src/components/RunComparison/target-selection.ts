// Pure state transitions of the charted-value selection: single pick, multi-value mode,
// and reconciliation when targets change. Transitions return the same object on a no-op.

import {CandidateOverrides, ComparisonTarget, candidateId} from './types';

export interface TargetSelection {
  selectedKey: string;
  multiMode: boolean;
  multiKeys: string[];
}

export const EMPTY_TARGET_SELECTION: TargetSelection = {selectedKey: '', multiMode: false, multiKeys: []};

export const isCharted = (state: TargetSelection, key: string) =>
  state.multiMode ? state.multiKeys.includes(key) : key === state.selectedKey;

/**
 * Row click. Single mode: selects, or clears when clicking the selection again; Shift+click
 * keeps a compatible selection and adds the clicked value, otherwise re-anchors on it and
 * enters multi mode only if it has partners. Multi mode: toggles the value; unchecking the
 * last one exits the mode with it selected, unchecking the anchor moves it to a charted value.
 */
export function clickTarget(
  state: TargetSelection,
  key: string,
  shift: boolean,
  compatibleKeys: (anchorKey: string) => string[],
): TargetSelection {
  const {selectedKey, multiMode, multiKeys} = state;
  if (!multiMode) {
    if (!shift)
      return {...state, selectedKey: selectedKey === key ? '' : key};
    if (selectedKey && compatibleKeys(selectedKey).includes(key))
      return {selectedKey, multiMode: true, multiKeys: [...new Set([selectedKey, key])]};
    return compatibleKeys(key).length > 1 ?
      {selectedKey: key, multiMode: true, multiKeys: [key]} :
      {...state, selectedKey: key};
  }
  const next = multiKeys.includes(key) ? multiKeys.filter((item) => item !== key) : [...multiKeys, key];
  if (next.length === 0)
    return {selectedKey: key, multiMode: false, multiKeys: []};
  return {selectedKey: next.includes(selectedKey) ? selectedKey : next[0], multiMode, multiKeys: next};
}

/** Mode toggle: entering seeds the selection with the current value, leaving keeps the first charted one. */
export function setMultiMode(state: TargetSelection, on: boolean): TargetSelection {
  if (on === state.multiMode)
    return state;
  if (on)
    return {...state, multiMode: true, multiKeys: state.selectedKey ? [state.selectedKey] : []};
  return {...state, multiMode: false, selectedKey: state.multiKeys[0] ?? state.selectedKey};
}

/**
 * Drops multi-value keys whose target vanished; merely conflicting targets stay selected.
 * Losing every key exits the mode; a vanished anchor moves to a still-charted value.
 */
export function reconcileTargets(state: TargetSelection, keys: Set<string>): TargetSelection {
  let next = state;
  const kept = state.multiKeys.filter((key) => keys.has(key));
  if (kept.length !== state.multiKeys.length)
    next = {...next, multiKeys: kept, multiMode: kept.length > 0 && next.multiMode};
  if (next.multiMode && !keys.has(next.selectedKey) && kept.length > 0)
    next = {...next, selectedKey: kept[0]};
  return next;
}

/** Drops overrides whose target or candidate no longer resolves. */
export function pruneOverrides(overrides: CandidateOverrides, targets: ComparisonTarget[]): CandidateOverrides {
  const validIds = new Map(targets.flatMap((target) => target.kind === 'column' ?
    [[target.key, new Set(target.candidates.map((candidate) => candidateId(candidate.binding)))] as const] : []));
  let changed = false;
  const next: CandidateOverrides = {};
  for (const [key, ids] of Object.entries(overrides)) {
    const valid = validIds.get(key);
    const kept = valid ? Object.entries(ids).filter(([id]) => valid.has(id)) : [];
    if (kept.length !== Object.keys(ids).length)
      changed = true;
    if (kept.length > 0)
      next[key] = Object.fromEntries(kept);
  }
  return changed ? next : overrides;
}
