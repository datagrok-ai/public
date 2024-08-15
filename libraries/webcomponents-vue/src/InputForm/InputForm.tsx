import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, isReactive, KeepAlive, PropType, ref, shallowRef, watch, watchEffect} from 'vue';
import type {InputFormT} from '@datagrok-libraries/webcomponents/src';
import {getValidators,
  injectInputBaseValidation,
  isInputBase,
  validate,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SYNC_FIELD} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {ValidationResult, Validator} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {FuncCallInputValidated, isFuncCallInputValidated} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-input-form': InputFormT
    }
  }
}

export const InputForm = defineComponent({
  name: 'InputForm',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
  },
  setup(props, {emit}) {
    let currentForm = undefined as undefined | DG.InputForm;

    let loadedValidators = shallowRef({});

    watch(() => props.funcCall, () => {
      console.log('props.funcCall changed', props.funcCall);

      getValidators(props.funcCall, SYNC_FIELD.INPUTS).then((res) => loadedValidators.value = res);
    })

    watch(loadedValidators, () => {
      runValidation();
    })

    const allParams = (funcCall: DG.FuncCall) =>
      [...funcCall.inputParams.values() ?? []].map((param) => param.name);

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);
      currentForm = event.detail;
      if (!currentForm) return;

      allParams(props.funcCall)
        .map((param) => currentForm!.getInput(param))
        .filter((input) => isInputBase(input))
        .forEach((input) => injectInputBaseValidation(input));
    };

    const inputChangedCb = (event: {detail: DG.EventData<DG.InputArgs>}) => {
      runValidation([event.detail.args.input.property.name]);
    };

    const runValidation = async (paramNames?: string[]) => {
      if (!currentForm) return;

      const paramsToValidate = paramNames ?? allParams(props.funcCall);

      const controller = new AbortController();

      const results = await validate({isRevalidation: false},
        paramsToValidate, controller.signal, SYNC_FIELD.INPUTS, {
          funcCall: props.funcCall,
        }, loadedValidators.value);

      Object.keys(results)
        .map((paramName) => currentForm!.getInput(paramName))
        .filter((input) => isFuncCallInputValidated(input))
        .forEach((input) => input.setValidation(results[input.property.name]));
    };

    return () => {
      const form = <dg-input-form
        funcCall={props.funcCall}
        onFormReplaced={formReplacedCb}
        onInputChanged={inputChangedCb}>
      </dg-input-form>;

      return form;
    };
  },
});
