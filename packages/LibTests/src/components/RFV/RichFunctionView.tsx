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
          
    const historyRef = ref(null as InstanceType<typeof History> | null);
    const helpRef = ref(null as InstanceType<typeof MarkDown> | null);
    const formRef = ref(null as HTMLElement | null);
    const dockRef = ref(null as InstanceType<typeof DockManager> | null);
    const handlePanelClose = (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;
    };

    const save = () => {
      if (!dockRef.value) return;

      dockRef.value.saveLayout();
    };

    const load = () => {
      if (!dockRef.value) return;

      dockRef.value.loadLayout();
    };

    return () => (
      <div class='w-full h-full flex'>
        <RibbonPanel>
          <IconFA
            name='pen'
            tooltip={formHidden.value ? 'Open inputs': 'Close inputs'}
            onClick={() => formHidden.value = !formHidden.value} 
            style={{'background-color': !formHidden.value ? 'var(--grey-1)': null}}
          />
          <IconFA
            name='play'
            tooltip='Run step'
            onClick={run} 
          />
          <IconFA
            name='save'
            tooltip='Save the layout'
            onClick={save}
          />
          <IconFA
            name='life-ring'
            tooltip='Load the layout'
            onClick={load}
          />
          { hasContextHelp.value && <IconFA 
            name='info' 
            tooltip={ helpHidden.value ? 'Open help panel' : 'Close help panel' }
            onClick={() => helpHidden.value = !helpHidden.value}
            style={{'background-color': !helpHidden.value ? 'var(--grey-1)': null}}
          /> }
          <IconFA 
            name='history' 
            tooltip='Open history panel' 
            onClick={() => historyHidden.value = !historyHidden.value}
            style={{'background-color': !historyHidden.value ? 'var(--grey-1)': null}}
          />
        </RibbonPanel>
        <DockManager onPanelClosed={handlePanelClose} ref={dockRef}>
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
              ref={historyRef}
            />: null }
          
          { !formHidden.value ?
            <div 
              class='flex flex-col p-2'
              dock-spawn-dock-type='left'
              dock-spawn-dock-ratio={0.2}
              title='Inputs'
              ref={formRef}
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
              ref={helpRef}
            /> : null 
          }
        </DockManager>
      </div>
    );
  },
});
