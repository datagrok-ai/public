import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, PropType, reactive, ref} from 'vue';
import type {InputFormT} from '@datagrok-libraries/webcomponents/src';
import {getValidators,
  injectInputBaseValidation,
  isInputBase,
  validate,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SYNC_FIELD} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {Validator} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {FuncCallInputValidated, isFuncCallInputValidated} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-input-form': InputFormT
    }
  }
}

const getFormSource = (form: DG.InputForm): DG.FuncCall => {
  return DG.toJs(form.source);
};

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

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);

      const form = event.detail;

      currentForm = form;

      if (!form) return;

      [...getFormSource(form).inputParams.values()].map((param) => param.name).forEach((param) => {
        const input = form.getInput(param);
        if (isInputBase(input))
          injectInputBaseValidation(input);
      });

      const loadedValidators = await getValidators(getFormSource(form), SYNC_FIELD.INPUTS);
      inputValidators.value = loadedValidators;

      runValidation();
    };

    const inputChangedCb = () => {
      runValidation();
    };

    const runValidation = async (paramNames?: string[]) => {
      if (!currentForm) return;

      const paramsToValidate = paramNames ?? allParams(getFormSource(currentForm));

      const controller = new AbortController();

      const results = await validate({isRevalidation: false},
        paramsToValidate, controller.signal, SYNC_FIELD.INPUTS, {
          funcCall: getFormSource(currentForm),
        }, inputValidators.value);

      Object.entries(results).forEach(([paramName, res]) => {
        const t = currentForm!.getInput(paramName);
        if (isFuncCallInputValidated(t))
          t.setValidation(res);
      });
    };

    const inputValidators = ref({} as Record<string, Validator>);

    getValidators(props.funcCall, SYNC_FIELD.INPUTS).then((loadedValidators) => {
      inputValidators.value = loadedValidators;

      runValidation();
    });

    return () => {
      const form = <dg-input-form
        funcCall={props.funcCall}
        onFormReplaced={formReplacedCb}
        onInputChanged={inputChangedCb}>
      </dg-input-form>;

      return (
        <keep-alive>
          { form }
        </keep-alive>
      );
    };
  },
});
