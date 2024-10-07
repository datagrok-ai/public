import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {RichFunctionView} from './RichFunctionView';
import {historyUtils} from '@datagrok-libraries/compute-utils';

export const RFVWrapper = Vue.defineComponent({
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
  },
  setup(props) {
    const currentFuncCall = Vue.shallowRef(props.funcCall);
    const currentCallState = Vue.ref(
      {isRunning: false, isOutputOutdated: false, isRunnable: true, pendingDependencies: []},
    );

    const func = Vue.computed(() => currentFuncCall.value.func);
    const isRunningOnInput = Vue.computed(() => func.value.options['runOnInput'] === 'true');

    const runFunc = async () => {
      if (!currentFuncCall.value) return;

      currentCallState.value.isRunning = true;
      await currentFuncCall.value.call();
      currentCallState.value.isRunning = false;

      if (!isRunningOnInput.value)     
        await historyUtils.saveRun(currentFuncCall.value);
    };

    return () => (
      <RichFunctionView 
        funcCall={currentFuncCall.value}
        callState={currentCallState.value}
        onUpdate:funcCall={(chosenCall) => currentFuncCall.value = chosenCall}
        onRunClicked={runFunc}
      />
    );
  },
});
