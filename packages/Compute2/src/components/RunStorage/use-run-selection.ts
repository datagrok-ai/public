import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

export type RunClick = Pick<MouseEvent, 'ctrlKey' | 'metaKey' | 'detail'>;

/** The runs selected in a list: a click selects one run, a ctrl-click adds or removes it. A run the list stops
 * showing is unselected, so every action sees only listed runs. */
export function useRunSelection(listed: () => DG.FuncCall[]) {
  const selection = Vue.shallowRef<DG.FuncCall[]>([]);
  const isSelected = (run: DG.FuncCall) => selection.value.some((item) => item.id === run.id);

  // a double click opens the run only when neither of its clicks toggled the selection
  let doubleClickToggles = false;
  /** Applies the click; true when it is a double click that opens the run. */
  const click = (run: DG.FuncCall, e: RunClick) => {
    const toggles = e.ctrlKey || e.metaKey;
    if (e.detail <= 1)
      doubleClickToggles = toggles;
    const toggled = isSelected(run) ?
      selection.value.filter((item) => item.id !== run.id) : [...selection.value, run];
    selection.value = toggles ? toggled : [run];
    return e.detail === 2 && !toggles && !doubleClickToggles;
  };

  /** The runs an action on the shown run applies to: all selected ones when it is among several selected. */
  const actionRuns = (shown: DG.FuncCall) =>
    selection.value.length > 1 && isSelected(shown) ? selection.value : [shown];

  Vue.watch(listed, (runs) => {
    const isListed = (run: DG.FuncCall) => runs.some((item) => item.id === run.id);
    if (!selection.value.every(isListed))
      selection.value = selection.value.filter(isListed);
  });

  const clear = () => selection.value = [];
  const keepOnly = (runs: DG.FuncCall[]) => selection.value = runs;

  return {selection: selection as Readonly<typeof selection>, isSelected, click, actionRuns, clear, keepOnly};
}
