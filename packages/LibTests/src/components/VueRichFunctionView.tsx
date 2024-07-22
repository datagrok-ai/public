import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, shallowRef, triggerRef, defineExpose, nextTick} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton} from '@datagrok-libraries/webcomponents-vue/src';
import {GridStack} from 'gridstack';
import 'gridstack/dist/gridstack.min.css';
import '../styles/VueRichFunctionView.css';

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
    let grid = null; 
    
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

    const viewerTypes = ref(['Scatter plot', 'Histogram', 'Line chart']);

    const runSimulation = async () => {
      await resolvedFuncCall.value.call();
      triggerRef(resolvedFuncCall);
    };

    const restoreDefault = () => {
      resolvedFuncCall.value = getDefaultCall();
    };

    let inited = false;

    onMounted(async () => {
      await nextTick();
      if (inited) return;

      const formNode = document.querySelector('#formNode')!;
        
      grid = GridStack.init({
        float: true,
        auto: false,
        column: 30,
        margin: 0,
      });
    
      grid.addWidget(formNode, {
        minW: Math.ceil(160 / grid.cellWidth()), 
        maxW: Math.ceil(formNode.getBoundingClientRect().width / grid.cellWidth()), 
        h: Math.ceil(formNode.getBoundingClientRect().height / grid.getCellHeight()),
      });
      for (let i = 1; i <= 3; i++) 
        grid.addWidget(document.querySelector(`#viewerNode${i}`)!, {h: 5, w: 12}); 
      
      inited = true;     
    });

    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <div class="grid-stack"></div>

        <div id="formNode">
          <InputForm funcCall={resolvedFuncCall.value}> </InputForm>
          <div style={{display: 'flex', position: 'sticky', bottom: '10px'}}>
            <BigButton onClick={restoreDefault}> Restore default </BigButton>
            <BigButton onClick={runSimulation}> Run </BigButton>
          </div>
        </div>
        <div id="viewerNode1">
          <div>
            <label for="viewer-choice1">Choose a viewer:</label>
            <select v-model={viewerTypes.value[0]} name="viewers" id="viewer-choice1">
              <option value="Scatter plot">Scatter plot</option>
              <option value="Line chart">Line chart</option>
              <option value="Grid">Grid</option>
              <option value="Histogram">Histogram</option>
            </select>
          </div>
          <Viewer 
            style={{width: '100%'}} 
            type={viewerTypes.value[0]}
            dataFrame={resolvedFuncCall.value.outputs['simulation']}> 
          </Viewer>
        </div>
        <div id="viewerNode2">
          <div>
            <label for="viewer-choice2">Choose a viewer:</label>
            <select v-model={viewerTypes.value[1]} name="viewers" id="viewer-choice2">
              <option value="Scatter plot">Scatter plot</option>
              <option value="Line chart">Line chart</option>
              <option value="Grid">Grid</option>
              <option value="Histogram">Histogram</option>
            </select>
          </div>
          <Viewer 
            style={{width: '100%'}} 
            type={viewerTypes.value[1]}
            dataFrame={resolvedFuncCall.value.outputs['simulation']}> 
          </Viewer>
        </div>
        <div id="viewerNode3">
          <div>
            <label for="viewer-choice3">Choose a viewer:</label>
            <select v-model={viewerTypes.value[2]} id="viewer-choice3" name="viewers">
              <option value="Scatter plot">Scatter plot</option>
              <option value="Line chart">Line chart</option>
              <option value="Grid">Grid</option>
              <option value="Histogram">Histogram</option>
            </select>
          </div>
          <Viewer 
            style={{width: '100%'}} 
            type={viewerTypes.value[2]}
            dataFrame={resolvedFuncCall.value.outputs['simulation']}> 
          </Viewer>
        </div>
      </div>
    );
  },
});
