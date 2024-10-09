import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {InputFormT} from '@datagrok-libraries/webcomponents';
import {getValidators,
  injectInputBaseValidation,
  isInputBase,
  validate,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SYNC_FIELD} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {isFuncCallInputValidated} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';
import { computedAsync } from '@vueuse/core'

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-input-form': InputFormT
    }
  }
}

export const InputForm = Vue.defineComponent({
  name: 'InputForm',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    validationEnabled: {
      type: Boolean,
      default: false,
    }
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
    'update:funcCall': (call: DG.FuncCall) => call,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);

    let currentForm = undefined as undefined | DG.InputForm;

    let loadedValidators = computedAsync(() => getValidators(currentCall.value, SYNC_FIELD.INPUTS), {}, {shallow: true});

    if (props.validationEnabled) {
        Vue.watch(loadedValidators, () => {
        runValidation();
      });
    }

    const allParams = (funcCall: DG.FuncCall) =>
      [...funcCall.inputParams.values() ?? []].map((param) => param.name);

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);
      currentForm = event.detail;
      if (!currentForm) return;

      allParams(currentCall.value)
        .map((param) => currentForm!.getInput(param))
        .filter((input) => isInputBase(input))
        .forEach((input) => injectInputBaseValidation(input));
    };

    const inputChangedCb = (event: {detail: DG.EventData<DG.InputArgs>}) => {
      if (props.validationEnabled) runValidation([event.detail.args.input.property.name]);
      emit('update:funcCall', currentCall.value)
    };

    const runValidation = async (paramNames?: string[]) => {
      if (!currentForm) return;

      const paramsToValidate = paramNames ?? allParams(currentCall.value);

      const controller = new AbortController();

      const results = await validate({isRevalidation: false},
        paramsToValidate, controller.signal, SYNC_FIELD.INPUTS, {
          funcCall: currentCall.value,
        }, loadedValidators.value);

      Object.keys(results)
        .map((paramName) => currentForm!.getInput(paramName))
        .filter((input) => isFuncCallInputValidated(input))
        .forEach((input) => input.setValidation(results[input.property.name]));
    };

    return () => <dg-input-form
      funcCall={currentCall.value}
      onFormReplaced={formReplacedCb}
      onInputChanged={inputChangedCb}>
    </dg-input-form>;
  },
});
