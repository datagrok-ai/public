import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {Advice, ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import type {InputFormT} from '@datagrok-libraries/webcomponents';
import {injectInputBaseValidation, isInputBase} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ValidationResultBase} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {
  FuncCallInputValidated,
  isFuncCallInputValidated,
} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';

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
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    },
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
    actionRequested: (actionUuid: string) => actionUuid,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);
    const validationStates = Vue.computed(() => props.validationStates);

    const convertNewValidationToOld = (newResult: ValidationResult): ValidationResultBase => {
      const convert = (error: Advice) => ({
        description: error.description,
        actions: error.actions?.map((actionItem) => ({actionName: actionItem.actionName, action: () => {
          emit('actionRequested', actionItem.action);
        }})),
      });
      return {
        errors: newResult.errors?.map(convert),
        warnings: newResult.warnings?.map(convert),
        notifications: newResult.notifications?.map(convert),
      };
    };

    let currentForm = undefined as undefined | DG.InputForm;

    const runValidations = () => {
      const newStates = validationStates.value;
      if (!newStates || !currentForm) return;

      [...currentCall.value.inputParams.values()]
        .filter((param) => (currentForm!.getInput(param.property.name)))
        .forEach((param) => {
          const input = currentForm!.getInput(param.property.name);
          const validationState = newStates[param.property.name];

          if (isFuncCallInputValidated(input)) {
            (input as FuncCallInputValidated).setValidation(
              validationState ? convertNewValidationToOld(validationState): undefined,
            );
          }
        });
    };

    Vue.watch(validationStates, () => {
      runValidations();
    });

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);
      currentForm = event.detail;

      if (!currentForm || !validationStates.value) return;

      [...currentCall.value.inputParams.values()]
        .filter((param) => isInputBase(currentForm!.getInput(param.property.name)))
        .forEach((param) => {
          const input = currentForm!.getInput(param.property.name);
          injectInputBaseValidation(input);
        });
      runValidations();
    };

    return () => <dg-input-form
      funcCall={currentCall.value}
      onFormReplaced={formReplacedCb}
      validationStates={validationStates.value}>
    </dg-input-form>;
  },
});
