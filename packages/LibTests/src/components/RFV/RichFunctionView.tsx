import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton} from '@datagrok-libraries/webcomponents-vue/src';
import {GridStack} from 'gridstack';
import 'gridstack/dist/gridstack.min.css';
import './RichFunctionView.css';
import {getPropViewers} from '@datagrok-libraries/compute-utils/shared-utils/utils';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const RichFunctionView = defineComponent({
  name: 'RichFunctionView',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall | string>,
      required: true,
    },
  },
  setup(props) {
    let grid = null; 
    
    const currentCall = computed(() => props.funcCall instanceof DG.FuncCall ?
      props.funcCall : DG.Func.byName(props.funcCall).prepare({
        ambTemp: 22,
        initTemp: 100,
        desiredTemp: 30,
        area: 0.06, 
        heatCap: 4200,
        heatTransferCoeff: 8.3,
        simTime: 21600,
        previousRun: null,
      }));

    const paramsWithViewers = computed(() => {
      return [
        ...currentCall.value.func.inputs,
        ...currentCall.value.func.outputs,
      ].filter((prop) =>
        prop.propertyType === DG.TYPE.DATA_FRAME && getPropViewers(prop).config.length !== 0,
      );
    });

    const viewers = ref([] as Element[]);
    const formNode = ref(null as null | Element);

    const addViewer = (el: Element) => viewers.value.push(el); 

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

    const runSimulation = async () => {
      await currentCall.value.call();
      triggerRef(currentCall);
    };

    let inited = false;

    onMounted(async () => {
      await nextTick();
      if (inited) return;
  
      grid = GridStack.init({
        float: true,
        auto: false,
        column: 30,
        margin: 0,
      });
    
      if (formNode.value) {
        grid.addWidget(formNode.value, {
          minW: Math.ceil(160 / grid.cellWidth()), 
          maxW: Math.ceil(formNode.value.getBoundingClientRect().width / grid.cellWidth()), 
          h: Math.ceil(formNode.value.getBoundingClientRect().height / grid.getCellHeight()),
        });
      }
      
      let idx = 0;
      while (document.querySelector(`#viewer${idx}`)) {
        grid.addWidget(document.querySelector(`#viewer${idx}`)!, {h: 5, w: 12});    
        idx++;
      }   
      
      inited = true;     
    });
          
    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <div class="grid-stack"></div>

        <div ref={formNode}>
          <InputForm funcCall={currentCall.value}> </InputForm>
          <div style={{display: 'flex', position: 'sticky', bottom: '10px'}}>
            <BigButton onClick={runSimulation}> Run </BigButton>
          </div>
        </div>
        { 
          paramsWithViewers.value
            .map((viewer) => getPropViewers(viewer))
            .map(({name, config: allConfigs}) => 
              allConfigs.map((options, idx) => 
                <div id={`viewer${idx.toString()}`}>
                  <Viewer
                    type={options['type'] as string}
                    options={options}
                    dataFrame={currentCall.value.inputs[name] ?? currentCall.value.outputs[name]}
                    style={{width: '100%'}} 
                  /> 
                </div>,
              ),
            )
        }
      </div>
    );
  },
});
