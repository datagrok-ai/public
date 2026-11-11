import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {DfExportEntry} from '../../project-export';
import {RenderStateItem, getEmptyTabToProperties, tabToProperties} from './rfv-tabs';

export function useRfvTabs(currentCall: Readonly<Vue.Ref<DG.FuncCall>>, deps: {
  callState: () => FuncCallStateInfo | undefined,
  hiddenByMeta: Readonly<Vue.Ref<Set<string>>>,
  formAsTab: Readonly<Vue.Ref<boolean>>,
  isOutputOutdated: Readonly<Vue.Ref<boolean | undefined>>,
}) {
  const {hiddenByMeta, formAsTab, isOutputOutdated} = deps;
  const tabsData = Vue.shallowRef<RenderStateItem[]>([]);
  const tabLabels = Vue.shallowRef<string[]>([]);
  const tabToPropertiesMap = Vue.shallowRef(getEmptyTabToProperties());
  const userClosed = Vue.shallowRef(new Set<string>());

  // Mounted viewer per tab label; entries of closed tabs go stale (no undefined emit on
  // unmount) but are never read — collectDfExportEntries iterates visibleTabLabels only
  const liveViewers = new Map<string, DG.Viewer>();

  const visibleTabLabels = Vue.computed(() => tabLabels.value.filter((label) => {
    if (userClosed.value.has(label))
      return false;
    const inputContent = tabToPropertiesMap.value.inputs.get(label);
    if (inputContent && inputContent.type === 'dataframe' && hiddenByMeta.value.has(inputContent.name))
      return false;
    return true;
  }));

  // tabToProperties subscribes to param changes via useObservable, which needs an
  // active effect scope to register disposal — watcher callbacks have none, so each
  // rebuild runs in its own scope and stopping it releases the previous build's
  // subscriptions (down to the Dart-side param listeners)
  let tabsScope: Vue.EffectScope | undefined;
  const rebuildTabs = (call: DG.FuncCall) => {
    tabsScope?.stop();
    tabsScope = Vue.effectScope();
    tabToPropertiesMap.value = tabsScope.run(() => tabToProperties(call))!;
    tabLabels.value = [
      ...tabToPropertiesMap.value.inputs.keys(),
      ...tabToPropertiesMap.value.outputs.keys(),
    ];
  };
  Vue.onBeforeUnmount(() => tabsScope?.stop());

  // Per-function preferred tab, tracked separately for the input and output sides and pushed
  // to the dock as `preferredPanelTitle`. Restored only on function switch and run completion
  // (see the watcher below); a user click sets it directly so the click is respected until
  // the next restore.
  const inputKey = () => `opened_input_tab_${currentCall.value?.func?.nqName}`;
  const outputKey = () => `opened_output_tab_${currentCall.value?.func?.nqName}`;
  const preferredTab = Vue.ref<string | null>(null);

  const sideTabs = (side: 'inputs' | 'outputs') =>
    visibleTabLabels.value.filter((l) => tabToPropertiesMap.value[side].has(l));

  // Saved tab if still visible, else the default for the side (not persisted): the last output
  // tab (typically the final result) but the first input tab.
  const resolveSide = (side: 'inputs' | 'outputs', key: string) => {
    const tabs = sideTabs(side);
    const saved = sessionStorage.getItem(key);
    if (saved && tabs.includes(saved))
      return saved;
    return (side === 'outputs' ? tabs[tabs.length - 1] : tabs[0]) ?? null;
  };

  // formAsTab forces the 'Inputs' tab; otherwise input vs output side by run state.
  const resolvePreferredTab = () =>
    formAsTab.value ? 'Inputs' :
      (isOutputOutdated.value ?
        resolveSide('inputs', inputKey()) :
        resolveSide('outputs', outputKey()));

  const handleTabClicked = (title: string | null) => {
    if (!title)
      return;
    preferredTab.value = title;
    // 'Inputs' (form tab or side-panel) is not persisted: forced by formAsTab or a sticky panel.
    if (tabToPropertiesMap.value.inputs.has(title))
      sessionStorage.setItem(inputKey(), title);
    else if (tabToPropertiesMap.value.outputs.has(title))
      sessionStorage.setItem(outputKey(), title);
  };

  // registered before the tab data watcher, which reads the rebuilt tabs and the reset closed set
  Vue.watch(currentCall, (call) => {
    rebuildTabs(call);
    userClosed.value = new Set();
  }, {immediate: true});

  Vue.watch([currentCall, deps.callState, visibleTabLabels], (
    [call, callState, labels], [prevCall, prevCallState, prevLabels],
  ) => {
    if (prevCall === call && prevCallState === callState && prevLabels === labels)
      return;
    // Re-run mutates outputs in place; refresh tabToPropertiesMap and tabLabels so
    // hide-empty-outputs can drop/re-add tabs whose emptiness changed.
    if (prevCall && call === prevCall && callState !== prevCallState)
      rebuildTabs(call);
    const map = tabToPropertiesMap.value;

    tabsData.value = labels.map((tabLabel) =>
      ({
        tabLabel,
        tabContent: map.inputs.get(tabLabel) ?? map.outputs.get(tabLabel)!,
        isInput: !!map.inputs.has(tabLabel),
      }));

    // Restore the preferred tab on function switch (incl. initial mount) and on run completion
    // (isOutputOutdated true->false). A plain visibleTabLabels change touches neither, so a
    // mid-step tab show/hide leaves focus untouched.
    const switched = prevCall !== call;
    const justRan = !switched && !!prevCallState?.isOutputOutdated && !!callState && !callState.isOutputOutdated;
    if (switched || justRan)
      preferredTab.value = resolvePreferredTab();
  }, {immediate: true});

  const closeTab = (label: string) => {
    if (!tabLabels.value.includes(label))
      return;
    const next = new Set(userClosed.value);
    next.add(label);
    userClosed.value = next;
  };

  const showAllTabs = () => {
    userClosed.value = new Set();
  };

  const trackViewer = (label: string, viewer: DG.Viewer | undefined) => {
    if (viewer) liveViewers.set(label, viewer);
    else liveViewers.delete(label);
  };

  // One entry per visible dataframe param, viewer configs merged across its tabs
  const collectDfExportEntries = (): DfExportEntry[] => {
    const map = tabToPropertiesMap.value;
    const entries = new Map<string, DfExportEntry>();
    for (const label of visibleTabLabels.value) {
      const isInput = map.inputs.has(label);
      const content = (isInput ? map.inputs : map.outputs).get(label);
      if (!content || content.type !== 'dataframe' || !content.df.value)
        continue;
      const key = `${isInput ? 'in' : 'out'}:${content.name}`;
      let entry = entries.get(key);
      if (!entry) {
        entry = {name: content.name, isInput, df: content.df.value, viewers: []};
        entries.set(key, entry);
      }
      const {type, ...annotated} = content.config;
      const live = liveViewers.get(label);
      const options = live ? {...annotated, ...live.getOptions().look} : annotated;
      entry.viewers.push({type: (type as string) ?? DG.VIEWER.GRID, options});
    }
    return [...entries.values()];
  };

  return {
    tabsData: Vue.shallowReadonly(tabsData),
    preferredTab: Vue.shallowReadonly(preferredTab),
    userClosed: Vue.shallowReadonly(userClosed),
    handleTabClicked, closeTab, showAllTabs, trackViewer, collectDfExportEntries,
  };
}
