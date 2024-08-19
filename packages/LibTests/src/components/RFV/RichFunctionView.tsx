import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button, TabHeaderStripe, Tabs, IconFA} from '@datagrok-libraries/webcomponents-vue/src';
import 'gridstack/dist/gridstack.min.css';
import './RichFunctionView.css';
import {categoryToDfParamMap, getPropViewers} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getDefaultValue} from '@datagrok-libraries/compute-utils/function-views/src/shared/utils';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const ScalarTable = defineComponent({
  name: 'ScalarTable',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
    category: {
      type: String,
      required: true,
    },
  },
  setup(props) {
    const categoryScalars = computed(() => {
      return [
        ...props.funcCall.func.inputs,
        ...props.funcCall.func.outputs,
      ].filter((prop) => 
        prop.propertyType !== DG.TYPE.DATA_FRAME && 
        prop.propertyType !== DG.TYPE.GRAPHICS &&
        prop.category === props.category,
      );
    });

    return () => 
      <table>
        { 
          categoryScalars.value.map((prop) => {
            const precision = prop.options.precision;

            const scalarValue = precision && 
                prop.propertyType === DG.TYPE.FLOAT && props.funcCall.outputs[prop.name] ?
              props.funcCall.outputs[prop.name].toPrecision(precision):
              props.funcCall.outputs[prop.name];
            const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;
            
            return <tr> 
              <td>
                { prop.caption ?? prop.name }{units}
              </td> 
              <td>
                { scalarValue ?? '[No value]' }
              </td> 
            </tr>;
          }) 
        }
      </table>;
  },
});

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

    const categoryToDfParam = computed(() => categoryToDfParamMap(currentCall.value.func));

    const tabLabels = computed(() => {
      return [
        ...Object.keys(categoryToDfParam.value.inputs),
        ...Object.keys(categoryToDfParam.value.outputs),
      ];
    });

    const selectedIdx = ref(0);

    const runSimulation = async () => {
      await currentCall.value.call();
      triggerRef(currentCall);
    };
          
    return () => (
      <div style={{width: '100%', height: '100%', display: 'flex'}}>
        <div style={{padding: '8px'}}> 
          <InputForm funcCall={currentCall.value}/>
          <div style={{display: 'flex', position: 'sticky', bottom: '0px'}}>
            <BigButton onClick={runSimulation}> Run </BigButton>
          </div>
        </div>
        <Tabs 
          items={tabLabels.value.map((label) => ({label}))} 
          selected={selectedIdx.value} 
          onUpdate:selected={(v) => selectedIdx.value = v}
          style={{width: '100%'}}
        >
          {{
            default: () =>
              tabLabels.value.map((tabLabel) => categoryToDfParam.value.inputs[tabLabel] ?? 
                categoryToDfParam.value.outputs[tabLabel])            
                .flatMap((tabProps) => tabProps.map((prop) => getPropViewers(prop)))
                .map(({name, config: allConfigs}) => 
                  <div style={{display: 'flex', flexDirection: 'column'}}>
                    {
                      allConfigs.map((options, idx) => 
                        <div key={`${currentCall.value.id}_${idx}`} style={{height: '300px'}}>
                          <Viewer
                            type={options['type'] as string}
                            options={options}
                            dataFrame={currentCall.value.inputs[name] ?? currentCall.value.outputs[name]}
                            style={{width: '100%'}} 
                          /> 
                        </div>,
                      )
                    }
                    <ScalarTable 
                      funcCall={currentCall.value} 
                      category={tabLabels.value[selectedIdx.value]}
                    />
                  </div>, 
                ),
            stripePostfix: () => <IconFA name='info' tooltip='Open help panel'/>,
          }}
        </Tabs>
      </div>
    );
  },
});
