import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, isReactive, KeepAlive, MaybeRefOrGetter, onBeforeMount, onBeforeUnmount, onUpdated, PropType, reactive, Ref, ref, toValue, watch, watchEffect} from 'vue';
import type {InputFormT} from '@datagrok-libraries/webcomponents/src';
import {getValidators,
  injectInputBaseValidation,
  isInputBase,
  validate,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {SYNC_FIELD} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {isFuncCallInputValidated} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-input-form': InputFormT
    }
  }
}

const useValidators = async (funcCall: DG.FuncCall | MaybeRefOrGetter<DG.FuncCall>, isInput: SYNC_FIELD) => {
  return getValidators(toValue(funcCall), isInput);
};

const allParams = (funcCall: DG.FuncCall) =>
  [...funcCall.inputParams.values() ?? []].map((param) => param.name);

export const InputForm = defineComponent({
  name: 'InputForm',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
  },
  data() {
    return {
      validators: {},
      currentForm: undefined as undefined | DG.InputForm,
    };
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
  },
  watch: {
    funcCall: {
      async handler(val) {
        this.validators = await useValidators(this.funcCall, SYNC_FIELD.INPUTS);
        console.log('old style watch called', this.validators);
      },
    },
  },
  methods: {
    async formReplacedCb(event: {detail?: DG.InputForm}) {
      this.$emit('formReplaced', event.detail);
      this.currentForm = event.detail;
      if (!this.currentForm) return;

      allParams(this.funcCall)
        .map((param) => this.currentForm!.getInput(param))
        .filter((input) => isInputBase(input))
        .forEach((input) => injectInputBaseValidation(input));

      this.runValidation();
    },
    inputChangedCb(event: {detail: DG.EventData<DG.InputArgs>}) {
      this.runValidation([event.detail.args.input.property.name]);
    },
    async runValidation(paramNames?: string[]) {
      if (!this.currentForm) return;

      const paramsToValidate = paramNames ?? allParams(this.funcCall);

      const controller = new AbortController();

      const results = await validate({isRevalidation: false},
        paramsToValidate, controller.signal, SYNC_FIELD.INPUTS, {
          funcCall: this.funcCall,
        }, this.validators);

      Object.keys(results)
        .map((paramName) => this.currentForm!.getInput(paramName))
        .filter((input) => isFuncCallInputValidated(input))
        .forEach((input) => input.setValidation(results[input.property.name]));
    },
  },
  render() {
    const form = <dg-input-form
      funcCall={this.funcCall}
      onFormReplaced={this.formReplacedCb}
      onInputChanged={this.inputChangedCb}>
    </dg-input-form>;

    return form;
  },
});
