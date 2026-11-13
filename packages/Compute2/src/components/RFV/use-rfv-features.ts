import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';

interface DockSpawnConfigItem {
  'dock-spawn-dock-type'?: 'left' | 'right' | 'up' | 'down',
  'dock-spawn-dock-to'?: string,
  'dock-spawn-dock-ratio'?: number,
}

export function useRfvFeatures(currentCall: Readonly<Vue.Ref<DG.FuncCall>>, standalone: Readonly<Vue.Ref<boolean>>) {
  const func = () => currentCall.value.func;
  const features = Vue.computed(() => Utils.getFeatures(func()));
  return {
    isSAenabled: Vue.computed(() => Utils.getFeature(features.value, 'sens-analysis', false)),
    isFittingEnabled: Vue.computed(() => Utils.getFeature(features.value, 'fitting', false)),
    allowRerun: Vue.computed(() => Utils.getFeature(features.value, 'rerun', false)),
    runLabel: Vue.computed(() => Utils.getRunLabel(func()) ?? 'Run'),
    dockSpawnConfig: Vue.computed(() => Utils.getDockSpawnConfig(func()) as Record<string, DockSpawnConfigItem>),
    formAsTab: Vue.computed(() => Utils.getFormAsTab(func())),
    autorun: Vue.computed(() => standalone.value && Utils.isRunningOnInput(func())),
    inputToggles: Vue.computed(() => standalone.value ? [] : Utils.getInputToggles(func())),
  };
}
