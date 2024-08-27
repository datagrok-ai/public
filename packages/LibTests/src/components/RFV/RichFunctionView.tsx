import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch, shallowRef} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button, TabHeaderStripe, Tabs, IconFA} from '@datagrok-libraries/webcomponents-vue/src';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';

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
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
  },
  setup(props, {emit}) {
    const currentCall = computed(() => props.funcCall);

    const categoryToDfParam = computed(() => Utils.categoryToDfParamMap(currentCall.value.func));

    const tabLabels = computed(() => {
      return [
        ...Object.keys(categoryToDfParam.value.inputs),
        ...Object.keys(categoryToDfParam.value.outputs),
      ];
    });

    const selectedIdx = ref(0);

    const run = async () => {
      currentCall.value.newId();
      await currentCall.value.call();
      emit('update:funcCall', currentCall.value);
    };

    const formHidden = ref(false);
    const historyHidden = ref(true);

    const hasContextHelp = computed(() => Utils.hasContextHelp(currentCall.value.func));
          
    return () => (
      <div style={{width: '100%', height: '100%', display: 'flex'}}>
        <div style={{display: 'flex', flexDirection: 'column'}}>
          <div style={{display: 'flex', flexDirection: 'row', justifyContent: 'flex-end', padding: '0px 8px'}}>
            <IconFA 
              name={formHidden.value ? 'chevron-double-right': 'chevron-double-left'} 
              tooltip={formHidden.value ? 'Show input form': 'Hide input form'}
              onClick={() => formHidden.value = !formHidden.value}
            />
          </div>
          <div style={{padding: '8px', display: formHidden.value ? 'none': 'block'}}> 
            <InputForm funcCall={currentCall.value}/>
            <div style={{display: 'flex', position: 'sticky', bottom: '0px'}}>
              <BigButton onClick={run}> Run </BigButton>
            </div>
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
                .flatMap((tabProps) => tabProps.map((prop) => Utils.getPropViewers(prop)))
                .map(({name, config: allConfigs}) => 
                  <div style={{display: 'flex', flexDirection: 'column'}}>
                    {
                      allConfigs.map((options) => 
                        <Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[name] ?? currentCall.value.outputs[name]}
                          style={{width: '100%', height: '300px'}} 
                        />, 
                      )
                    }
                    <ScalarTable 
                      funcCall={currentCall.value} 
                      category={tabLabels.value[selectedIdx.value]}
                    />
                  </div>),
            stripePostfix: () => <div style={{display: 'flex', gap: '6px', padding: '0px 6px'}}>
              { hasContextHelp.value && <IconFA 
                name='info' 
                tooltip='Open help panel' 
                onClick={async () => Utils.showHelpWithDelay((await Utils.getContextHelp(currentCall.value.func))!)}
              /> }
              <IconFA 
                name='history' 
                tooltip='Open history panel' 
                onClick={() => historyHidden.value = !historyHidden.value}
              />
            </div>,
          }}
        </Tabs>
        <History 
          func={currentCall.value.func}
          showActions
          showBatchActions
          isHistory
          style={{display: historyHidden.value ? 'none': 'block', width: '30%'}}
          onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
        />
      </div>
    );
  },
});
