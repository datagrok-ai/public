import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {Advice, ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import type {InputFormT} from '@datagrok-libraries/webcomponents';
import {injectInputBaseValidation, injectLockStates, isInputBase} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ValidationResultBase} from '@datagrok-libraries/compute-utils/shared-utils/validation';
import {
  FuncCallInput,
  FuncCallInputValidated,
  isFuncCallInputValidated,
  isInputLockable,
} from '@datagrok-libraries/compute-utils/shared-utils/input-wrappers';
import $ from 'cash-dom';
import {injectLockIcons} from '@datagrok-libraries/compute-utils/function-views/src/shared/utils';
import {ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FuncCallInputStatusable, injectInputBaseStatus, isInputInjected} from './utils';
import {BehaviorSubject} from 'rxjs';

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
    consistencyStates: {
      type: Object as Vue.PropType<Record<string, ConsistencyInfo>>,
    },
    callMeta: {
      type: Object as Vue.PropType<Record<string, BehaviorSubject<any>>>,
    },
    isReadonly: {
      type: Boolean,
    },
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
    actionRequested: (actionUuid: string) => actionUuid,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);
    const validationStates = Vue.computed(() => props.validationStates);
    const consistencyStates = Vue.computed(() => props.consistencyStates);
    const isReadonly = Vue.computed(() => props.isReadonly);

    const currentForm = Vue.ref(undefined as undefined | DG.InputForm);

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

    const callMeta = Vue.computed(() => props.callMeta);

    Vue.watchEffect(() => {
      if (!currentForm.value) return;

      [...currentCall.value.inputParams.values()]
        .filter((param) => (currentForm.value!.getInput(param.property.name)))
        .forEach((param) => {
          const input = currentForm.value!.getInput(param.property.name);
          const validationState = validationStates.value?.[param.property.name];
          const consistencyState = consistencyStates.value?.[param.property.name];
          const paramItems = callMeta.value?.[param.property.name].value?.['items'];

          if (!isInputInjected(input))
            injectInputBaseStatus(input);

          (input as any).setStatus({
            validation: validationState ? convertNewValidationToOld(validationState): undefined,
            consistency: consistencyState,
          });

          input.enabled = !isReadonly.value;

          if (paramItems && input.inputType === DG.InputType.Choice)
            (input as DG.ChoiceInput<any>).items = paramItems;
        });
    });

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);
      currentForm.value = event.detail;
    };

    return () => <dg-input-form
      funcCall={currentCall.value}
      onFormReplaced={formReplacedCb}>
    </dg-input-form>;
  },
});
