import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, ref, shallowRef, triggerRef, watch} from 'vue';
import {RichFunctionView} from './RichFunctionView';
import {Button} from '@datagrok-libraries/webcomponents-vue';
import {deepCopy} from '@datagrok-libraries/compute-utils/shared-utils/utils';

export const RFVTestApp = defineComponent({
  setup() {
    const initialName = 'Compute:ObjectCooling';

    const currentFuncCall = shallowRef(DG.Func.byName(initialName).prepare());

    const changeFunc = () => {
      const nfc = currentFuncCall.value.func.nqName === 'LibTests:TestAdd2' ? 
        'Compute:ObjectCooling':
        'LibTests:TestAdd2';
      currentFuncCall.value = DG.Func.byName(nfc).prepare();
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <Button onClick={changeFunc}> Change func </Button>
        <RichFunctionView 
          funcCall={currentFuncCall.value}
          onUpdate:funcCall={(chosenCall) => currentFuncCall.value = chosenCall}
        />
      </div>
    );
  },
});
