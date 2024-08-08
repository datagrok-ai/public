import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch, shallowRef} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button} from '@datagrok-libraries/webcomponents-vue/src';
import {GridStack} from 'gridstack';
import wu from 'wu';
import 'gridstack/dist/gridstack.min.css';
import './RichFunctionView.css';
import {getPropViewers} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getDefaultValue} from '@datagrok-libraries/compute-utils/function-views/src/shared/utils';

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
    let grid = null as null | GridStack; 
    
    const currentCall = computed(() => {
      if (props.funcCall instanceof DG.FuncCall) 
        return props.funcCall;

      const func = DG.Func.byName(props.funcCall);
      return func.prepare(func.inputs.reduce((acc, prop) => {
        acc[prop.name] = getDefaultValue(prop);
        return acc;
      }, {} as Record<string, any>));
    });

    watch(currentCall, async (newVal, oldVal) => {
      if (newVal.func.id === oldVal.func.id) return;
      
      grid?.removeAll(true);
      await nextTick();
      populateGridStack();
    });

    const paramsWithViewers = computed(() => [
      ...currentCall.value.func.inputs,
      ...currentCall.value.func.outputs,
    ].filter((prop) =>
      prop.propertyType === DG.TYPE.DATA_FRAME && getPropViewers(prop).config.length !== 0,
    ));

    const formNode = ref(null as null | Element);

    const populateGridStack = () => {   
      if (!grid) return;
      
      if (formNode.value) {
        grid.addWidget(formNode.value, {
          minW: Math.ceil(160 / grid.cellWidth()), 
          maxW: Math.ceil(formNode.value.getBoundingClientRect().width / grid.cellWidth()), 
          h: Math.ceil(formNode.value.getBoundingClientRect().height / grid.getCellHeight()),
        });
      }

      document.querySelectorAll(`[id^='viewer']`).forEach(
        (node) => grid!.addWidget(node as HTMLElement, {h: 5, w: 12}),
      );
    };

    const runSimulation = async () => {
      await currentCall.value.call();
      triggerRef(currentCall);
    };

    onMounted(async () => {
      await nextTick();

      grid = GridStack.init({
        float: true,
        auto: false,
        column: 30,
        margin: 0,
      });
       
      populateGridStack();
    });
          
    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <div class="grid-stack" />
        <div ref={formNode}>
          <InputForm funcCall={currentCall.value} />
          <div style={{display: 'flex', position: 'sticky', bottom: '10px'}}>
            <Button> Reset </Button>
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
