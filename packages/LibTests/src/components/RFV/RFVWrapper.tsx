import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {RichFunctionView} from './RichFunctionView';

export const RFVWrapper = Vue.defineComponent({
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
  },
  setup(props) {
    const currentFuncCall = Vue.shallowRef(props.funcCall);
    const runFunc = async () => {
      if (!currentFuncCall.value) return;

      await currentFuncCall.value.call();
    };

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <RichFunctionView 
          funcCall={currentFuncCall.value}
          onUpdate:funcCall={(chosenCall) => currentFuncCall.value = chosenCall}
          onRunClicked={runFunc}
        />
      </div>
    );
  },
});
