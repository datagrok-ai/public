import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {InputFormT} from '@datagrok-libraries/webcomponents';
import {
  injectInputBaseValidation,
  isInputBase,
} from '@datagrok-libraries/compute-utils/shared-utils/utils';

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
  },
  emits: {
    formReplaced: (a: DG.InputForm | undefined) => a,
    'update:funcCall': (call: DG.FuncCall) => call,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);

    const inputChangedCb = async (event: {detail: DG.EventData<DG.InputArgs>}) => {
      currentCall.value.inputs[event.detail.args.input.property.name] = event.detail.args.input.value;
      emit('update:funcCall', currentCall.value)
    };

    return () => <dg-input-form
      funcCall={currentCall.value}
      onInputChanged={inputChangedCb}>
    </dg-input-form>;
  },
});
