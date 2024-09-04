import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch, shallowRef} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Tabs, IconFA, RibbonPanel, DockManager, MarkDown} from '@datagrok-libraries/webcomponents-vue/src';
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
        ...props.funcCall.func.inputs
          .filter((prop) => prop.category === props.category),
        ...props.funcCall.func.outputs
          .filter((prop) => 
            prop.category === props.category ||
            (['Misc', 'Output'].includes(prop.category) && props.category === 'Output'),
          ),
      ].filter((prop) => 
        prop.propertyType !== DG.TYPE.DATA_FRAME && 
        prop.propertyType !== DG.TYPE.GRAPHICS,
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
    const currentCall = computed(() => Utils.deepCopy(props.funcCall));

    const categoryToDfParam = computed(() => Utils.categoryToDfParamMap(currentCall.value.func));

    const tabLabels = computed(() => {
      return [
        ...Object.keys(categoryToDfParam.value.inputs),
        ...Object.keys(categoryToDfParam.value.outputs),
      ];
    });

    const run = async () => {
      currentCall.value.newId();
      await currentCall.value.call();
      emit('update:funcCall', currentCall.value);
    };

    const formHidden = ref(false);
    const historyHidden = ref(true);
    const helpHidden = ref(true);

    const hasContextHelp = computed(() => Utils.hasContextHelp(currentCall.value.func));

    const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';

    const helpText = ref(null as null | string);
    watch(currentCall, async () => {
      const loadedHelp = await Utils.getContextHelp(currentCall.value.func);

      helpText.value = loadedHelp ?? null;
    }, {immediate: true});
          
    return () => (
      <div class='w-full h-full flex'>
        <RibbonPanel>
          { hasContextHelp.value && <IconFA 
            name='info' 
            tooltip={ helpHidden.value ? 'Open help panel' : 'Close help panel' }
            onClick={() => helpHidden.value = !helpHidden.value}
          /> }
          <IconFA 
            name='history' 
            tooltip='Open history panel' 
            onClick={() => historyHidden.value = !historyHidden.value}
          />
        </RibbonPanel>
        <DockManager>
          { !historyHidden.value ? 
            <History 
              func={currentCall.value.func}
              showActions
              showBatchActions
              isHistory
              onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.2}
              {...{title: 'History'}}
            />: null }
          
          { !formHidden.value ?
            <div 
              class='flex flex-col p-2'
              dock-spawn-dock-type='left'
              dock-spawn-dock-ratio={0.2}
              title='Inputs'
            >
              <InputForm funcCall={currentCall.value}/>
              <div class='flex sticky bottom-0'>
                <BigButton onClick={run}> Run </BigButton>
              </div>
            </div>: null }
          
          { 
            tabLabels.value
              .map((tabLabel) => ({tabLabel, tabDfProps: categoryToDfParam.value.inputs[tabLabel] ?? 
                categoryToDfParam.value.outputs[tabLabel]}))            
              .map(({tabLabel, tabDfProps}) => {
                return <div 
                  class='flex flex-col'
                  title={tabLabel}
                >
                  { tabDfProps.map((tabProp) => {
                    const allConfigs = Utils.getPropViewers(tabProp).config;

                    return allConfigs.map((options) => (            
                      <div class='flex flex-col h-1/2 pl-2'>
                        <h2> { dfBlockTitle(tabProp) } </h2>
                        <Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[tabProp.name] ?? currentCall.value.outputs[tabProp.name]}
                          class='w-full'
                          style={{'height': '300px'}} 
                        />
                      </div>));
                  })
                  }
                  <ScalarTable 
                    funcCall={currentCall.value} 
                    category={tabLabel}
                  />
                </div>;
              })
          }
          { !helpHidden.value && helpText.value ? 
            <MarkDown 
              markdown={helpText.value}
              {...{title: 'Help'}}
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.15}
            /> : null 
          }
        </DockManager>
      </div>
    );
  },
});
