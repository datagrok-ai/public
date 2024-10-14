import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import type {InputFormT} from '@datagrok-libraries/webcomponents';
import {computedAsync} from '@vueuse/core'

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
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);

    let currentForm = undefined as undefined | DG.InputForm;

    const formReplacedCb = async (event: {detail?: DG.InputForm}) => {
      emit('formReplaced', event.detail);
      currentForm = event.detail;
      if (!currentForm) return;
    };

    return () => <dg-input-form
      funcCall={currentCall.value}
      onFormReplaced={formReplacedCb}>
    </dg-input-form>;
  },
});
