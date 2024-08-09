import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch, shallowRef, onActivated, onDeactivated} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button} from '@datagrok-libraries/webcomponents-vue/src';
import {GridItemHTMLElement, GridStack, GridStackElement, GridStackOptions, GridStackWidget} from 'gridstack';
import wu from 'wu';
import 'gridstack/dist/gridstack.min.css';
import './RichFunctionView.css';
import {getPropViewers} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getDefaultValue} from '@datagrok-libraries/compute-utils/function-views/src/shared/utils';
import {from} from '@vueuse/rxjs';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

const layoutCache = {} as Record<string, GridStackOptions>;

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
      
      if (!grid) return;

      layoutCache[oldVal.id] = grid.save(false, true) as GridStackOptions;
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

    const formNode = ref(null as null | HTMLElement);

    const onFormReplaced = () => {
      console.log('onFormReplaced');
      if (!grid || !formNode.value) return;
        
      if (layoutCache[currentCall.value.id]?.children?.[0]) return;

      const maxHeight = Math.ceil(
        formNode.value.firstElementChild!.getBoundingClientRect().height / grid.getCellHeight(),
      ) + 2;
      const maxWidth = Math.ceil(
        formNode.value.firstElementChild!.getBoundingClientRect().width / grid.cellWidth(),
      ) + 2;
      grid.update(formNode.value, {
        maxW: maxWidth,
        w: maxWidth, 
        maxH: maxHeight,
        h: maxHeight,
      });
    };

    const populateGridStack = () => {  
      if (!grid) return;
      
      if (formNode.value) {    
        grid.addWidget(formNode.value,
          layoutCache[currentCall.value.id]?.children?.[0] ??
          {
            w: 10,
            h: 10,
          });
      }

      document.querySelectorAll(`[id^='viewer']`).forEach(
        (node, idx) => grid!.addWidget(
          node as HTMLElement, 
          layoutCache[currentCall.value.id]?.children?.[idx + 1] ?? {h: 5, w: 12},
        ),
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
        margin: 8,
      });
       
      populateGridStack();
    });
          
    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <div class="grid-stack" />
        <div ref={formNode}>
          <InputForm funcCall={currentCall.value} onFormReplaced={onFormReplaced}/>
          <div style={{display: 'flex', position: 'sticky', bottom: '0px'}}>
            <Button> Reset </Button>
            <BigButton onClick={runSimulation}> Run </BigButton>
          </div>
        </div>
        { 
          paramsWithViewers.value
            .map((prop) => getPropViewers(prop))
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
