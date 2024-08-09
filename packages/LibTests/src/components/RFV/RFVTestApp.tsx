import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, ref, shallowRef, triggerRef} from 'vue';
import {RichFunctionView} from './RichFunctionView';
import {Button} from '@datagrok-libraries/webcomponents-vue/src';

export const RFVTestApp = defineComponent({
  setup() {
    const currentFuncCallName = ref('Compute:ObjectCooling');

    const currentFuncCall = computed(() => {
      return DG.Func.byName(currentFuncCallName.value).prepare();
    });

    const changeFunccall = () => {
      const nfc = currentFuncCallName.value === 'LibTests:SimpleInputs' ? 
        'Compute:ObjectCooling':
        'LibTests:SimpleInputs';
      currentFuncCallName.value = nfc;
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <Button onClick={changeFunccall}> Change funccall </Button>
        <RichFunctionView 
          funcCall={currentFuncCall.value}
        />
      </div>
    );
  },
});
