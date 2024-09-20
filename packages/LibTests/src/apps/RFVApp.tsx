import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';

import {RichFunctionView} from '../components/RFV/RichFunctionView';
import {Button} from '@datagrok-libraries/webcomponents-vue';

export const RFVApp = Vue.defineComponent({
  setup() {
    const initialName = 'Compute:ObjectCooling';

    const currentFuncCall = Vue.shallowRef(DG.Func.byName(initialName).prepare());

    const changeFunc = () => {
      const nfc = currentFuncCall.value.func.nqName === 'LibTests:TestAdd2' ? 
        'Compute:ObjectCooling':
        'LibTests:TestAdd2';
      currentFuncCall.value = DG.Func.byName(nfc).prepare();
    };

    const runFunc = async () => {
      await currentFuncCall.value.call();
      currentFuncCall.value = Utils.deepCopy(currentFuncCall.value);
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <Button onClick={changeFunc}> Change func </Button>
        <RichFunctionView 
          funcCall={currentFuncCall.value}
          onUpdate:funcCall={(chosenCall) => currentFuncCall.value = chosenCall}
          onRunClicked={runFunc}
        />
      </div>
    );
  },
});
