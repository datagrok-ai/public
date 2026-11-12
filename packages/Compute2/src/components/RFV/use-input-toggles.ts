import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {tryParseJson} from '@datagrok-libraries/utils/src/string-utils';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {getToggledInputs} from './input-toggles';

export function useInputToggles(
  currentCall: Readonly<Vue.Ref<DG.FuncCall>>,
  inputToggles: Readonly<Vue.Ref<Utils.InputToggle[]>>,
  hiddenByMeta: Readonly<Vue.Ref<Set<string>>>,
  validationState: Readonly<Vue.Ref<Record<string, ValidationResult> | undefined>>,
  consistencyState: Readonly<Vue.Ref<Record<string, ConsistencyInfo> | undefined>>,
) {
  const shownToggles = Vue.shallowRef(new Set<string>());
  const shownTogglesKey = () => `shown_input_toggles_${currentCall.value?.func?.nqName}`;
  const switchToggle = (caption: string) => {
    const shown = new Set(shownToggles.value);
    if (!shown.delete(caption))
      shown.add(caption);
    shownToggles.value = shown;
    sessionStorage.setItem(shownTogglesKey(), JSON.stringify([...shown]));
  };

  Vue.watch(currentCall, () => {
    shownToggles.value = new Set(tryParseJson(sessionStorage.getItem(shownTogglesKey()) ?? '[]') ?? []);
  }, {immediate: true});

  const toggledInputs = Vue.computed(() => getToggledInputs(
    [...currentCall.value.inputParams.values()].map((param) => param.property),
    inputToggles.value,
    shownToggles.value,
    hiddenByMeta.value,
    validationState.value,
    consistencyState.value,
  ));

  return {shownToggles: Vue.shallowReadonly(shownToggles), switchToggle, toggledInputs};
}
