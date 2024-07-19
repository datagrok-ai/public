import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, PropType, ref, shallowRef, triggerRef} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, SplitH} from '@datagrok-libraries/webcomponents-vue/src';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const VueRichFunctionView = defineComponent({
  name: 'VueRichFunctionView',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall | string>,
      required: true,
    },
  },
  setup(props) {
    const getDefaultCall = () => props.funcCall instanceof DG.FuncCall ?
      props.funcCall : DG.Func.byName(props.funcCall).prepare({
        ambTemp: 22,
        initTemp: 100,
        desiredTemp: 30,
        area: 0.06, 
        heatCap: 4200,
        heatTransferCoeff: 8.3,
        simTime: 21600,
        previousRun: null,
      });


    const resolvedFuncCall = shallowRef(getDefaultCall());

    const viewerType = ref('Scatter plot');

    const runSimulation = async () => {
      await resolvedFuncCall.value.call();
      triggerRef(resolvedFuncCall);
    };

    const restoreDefault = () => {
      resolvedFuncCall.value = getDefaultCall();
    };

    return () => (
      <SplitH resize={true}>
        <div>
          <InputForm funcCall={resolvedFuncCall.value}> </InputForm>
          <div style={{display: 'flex'}}>
            <BigButton onClick={restoreDefault}> Restore default </BigButton>
            <BigButton onClick={runSimulation}> Run </BigButton>
          </div>
        </div>
        <div style={{width: '100%', paddingLeft: '10px'}}>
          <label for="viewer-choice">Choose a viewer:</label>
          <select v-model={viewerType.value} id="viewer-choice" name="viewers">
            <option value="Scatter plot">Scatter plot</option>
            <option value="Line chart">Line chart</option>
            <option value="Grid">Grid</option>
            <option value="Histogram">Histogram</option>
          </select>
          <Viewer 
            style={{height: '100%'}} 
            type={viewerType.value}
            dataFrame={resolvedFuncCall.value.outputs['simulation']}> 
          </Viewer>
        </div>
      </SplitH>
    );
  },
});
