import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent} from 'vue';
import type {FormT} from '@datagrok-libraries/webcomponents';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-form': FormT
    }
  }
}

export const Form = defineComponent({
  name: 'Form',
  props: {
    funcCall: DG.FuncCall,
  },
  emits: {
    formChanged: (a: DG.Form) => a,
  },
  setup(props, {emit}) {
    const formChangedCb = (event: any) => {
      emit('formChanged', event.detail);
    };
    return () => {
      const form = <dg-form funcCall={props.funcCall} onFormChanged={formChangedCb}></dg-form>;
      return (
        <keep-alive>
          { form }
        </keep-alive>
      );
    };
  },
});
