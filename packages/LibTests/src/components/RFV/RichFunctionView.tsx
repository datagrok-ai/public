import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button} from '@datagrok-libraries/webcomponents-vue/src';
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
    const currentCall = computed(() => {
      if (props.funcCall instanceof DG.FuncCall) 
        return props.funcCall;

      const func = DG.Func.byName(props.funcCall);
      return func.prepare(func.inputs.reduce((acc, prop) => {
        acc[prop.name] = getDefaultValue(prop);
        return acc;
      }, {} as Record<string, any>));
    });

    const paramsWithViewers = computed(() => [
      ...currentCall.value.func.inputs,
      ...currentCall.value.func.outputs,
    ].filter((prop) =>
      prop.propertyType === DG.TYPE.DATA_FRAME && getPropViewers(prop).config.length !== 0,
    ));

    const runSimulation = async () => {
      await currentCall.value.call();
      triggerRef(currentCall);
    };
          
    return () => (
      <div style={{width: '100%', height: '100%'}}>
        <div>
          <InputForm funcCall={currentCall.value}/>
          <div style={{display: 'flex', position: 'sticky', bottom: '0px'}}>
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
