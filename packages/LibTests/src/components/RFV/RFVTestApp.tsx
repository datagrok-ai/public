import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, ref, shallowRef, triggerRef, watch} from 'vue';
import {RichFunctionView} from './RichFunctionView';
import {Button} from '@datagrok-libraries/webcomponents-vue/src';

export const RFVTestApp = defineComponent({
  setup() {
    const initialName = 'Compute:ObjectCooling';

    const currentFuncCall = shallowRef(DG.Func.byName(initialName).prepare());

    const changeFunc = () => {
      const nfc = currentFuncCall.value.func.name === 'LibTests:SimpleInputs' ? 
        'Compute:ObjectCooling':
        'LibTests:SimpleInputs';
      currentFuncCall.value.func.name = nfc;
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <Button onClick={changeFunc}> Change funccall </Button>
        <RichFunctionView 
          funcCall={currentFuncCall.value}
          onFuncCallChange={(chosenCall) => currentFuncCall.value = chosenCall}
        />
      </div>
    );
  },
});
