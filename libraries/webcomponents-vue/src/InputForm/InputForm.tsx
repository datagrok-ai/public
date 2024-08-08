import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, isReactive, KeepAlive, PropType, ref, watch, watchEffect} from 'vue';
import type {InputFormT} from '@datagrok-libraries/webcomponents/src';
import {getValidators,
  injectInputBaseValidation,
  isInputBase,
  validate,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SYNC_FIELD} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {Validator} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {isFuncCallInputValidated} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';
import {computedAsync} from '@vueuse/core';

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

    const allParams = (funcCall: DG.FuncCall) =>
      [...funcCall.inputParams.values() ?? []].map((param) => param.name);

    let loadedValidators = {};

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);

      const form = event.detail;

      currentForm = form;

      if (!form) return;

      [...props.funcCall.inputParams.values()]
        .map((param) => param.name)
        .forEach((param) => {
          const input = form.getInput(param);
          if (isInputBase(input))
            injectInputBaseValidation(input);
        });

      runValidation();
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
        }, loadedValidators);

      Object.entries(results).forEach(([paramName, res]) => {
        const t = currentForm!.getInput(paramName);
        if (isFuncCallInputValidated(t))
          t.setValidation(res);
      });
    };

    return () => {
      getValidators(props.funcCall, SYNC_FIELD.INPUTS).then((res) => loadedValidators = res);

      const form = <dg-input-form
        funcCall={props.funcCall}
        onFormReplaced={formReplacedCb}
        onInputChanged={inputChangedCb}>
      </dg-input-form>;

      return form;
    };
  },
});
