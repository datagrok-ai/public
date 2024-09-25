import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {type ViewerT} from '@datagrok-libraries/webcomponents';
import {
  Viewer, InputForm, 
  BigButton, IconFA, 
  RibbonPanel, DockManager, MarkDown,
  ComboPopup,
  RibbonMenu,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {useStorage} from '@vueuse/core';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';

type PanelsState = {
  historyHidden: boolean,
  helpHidden: boolean,
  formHidden: boolean,
  visibleTabLabels: string[],
};

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const ScalarsPanel = Vue.defineComponent({
  name: 'ScalarsPanel',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    category: {
      type: String,
      required: true,
    },
  },
  setup(props) {
    const categoryScalars = Vue.computed(() => {
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

    const getContent = (prop: DG.Property) => {
      const precision = prop.options.precision;

      const scalarValue = precision && 
                prop.propertyType === DG.TYPE.FLOAT && props.funcCall.outputs[prop.name] ?
        props.funcCall.outputs[prop.name].toPrecision(precision):
        props.funcCall.outputs[prop.name];
      const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

      return [scalarValue, units];
    };

    return () => 
      categoryScalars.value.length <= 3 ?
        <div 
          class='flex flex-wrap justify-around' 
          dock-spawn-dock-type='down'
          dock-spawn-dock-ratio={0.15}
        >
          { categoryScalars.value.map((prop) => {
            const [scalarValue, units] = getContent(prop);

            return <div 
              class='flex flex-col p-2 items-center gap-4 flex-nowrap'
            > 
              <div class='text-center' style={{color: 'var(--grey-4)'}}> { prop.caption ?? prop.name } </div>
              <span style={{fontSize: 'var(--font-size-large)'}}> {scalarValue} {units} </span>
            </div>;
          })}
        </div> :
        <div 
          class='flex flex-col border-t-1 border-[#f2f2f5]' 
          dock-spawn-dock-type='fill'
        >
          { 
            categoryScalars.value.map((prop) => {
              const [scalarValue, units] = getContent(prop);
            
              return <div class='
                flex justify-between items-center h-8 px-2 
                hover:bg-[#f2f2f5]'
              style={{borderBottom: 'solid #f2f2f5'}}
              > 
                <div class='flex flex-col'>
                  { prop.caption ?? prop.name }
                  <span style={{color: 'var(--grey-4)'}}> { prop.description && prop.description } </span>
                </div> 
                <div>
                  <span> 
                    { scalarValue ?? '[No value]' } {units} 
                  </span>
                </div> 
              </div>;
            }) 
          }
        </div>;
  },
});

export const RichFunctionView = Vue.defineComponent({
  name: 'RichFunctionView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    callState: {
      type: Object as Vue.PropType<FuncCallStateInfo>,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
    'runClicked': () => {},
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => props.funcCall);

    const categoryToDfParam = Vue.computed(() => Utils.categoryToDfParamMap(currentCall.value.func));

    const tabLabels = Vue.computed(() => {
      return [
        ...Object.keys(categoryToDfParam.value.inputs),
        ...Object.keys(categoryToDfParam.value.outputs),
      ];
    });

    const run = async () => {
      emit('runClicked');
    };

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);

    const hasContextHelp = Vue.computed(() => Utils.hasContextHelp(currentCall.value.func));

    const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';

    const helpText = Vue.ref(null as null | string);
    Vue.watch(currentCall.value, async () => {
      const loadedHelp = await Utils.getContextHelp(currentCall.value.func);

      helpText.value = loadedHelp ?? null;
    }, {immediate: true});
    
    const root = Vue.ref(null as HTMLElement | null);
    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);
    const formRef = Vue.shallowRef(null as HTMLElement | null);
    const dockRef = Vue.shallowRef(null as InstanceType<typeof DockManager> | null);
    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const tabIdx = visibleTabLabels.value.findIndex((label) => label === el.title);
      if (tabIdx >= 0) 
        visibleTabLabels.value.splice(tabIdx, 1);  
    };

    const saveLayout = () => {
      if (!dockRef.value) return;

      panelsState.value = JSON.stringify({
        historyHidden: historyHidden.value,
        helpHidden: helpHidden.value,
        formHidden: formHidden.value,
        visibleTabLabels: visibleTabLabels.value,
      });
      dockRef.value.saveLayout();
    };

    const panelsStorageName = Vue.computed(() => `${currentCall.value.func.nqName}_panels`);

    const panelsState = useStorage(
      panelsStorageName.value,
      null as null | string,
    );

    const loadLayout = async () => {
      if (!dockRef.value || !panelsState.value) return;

      const openedPanels = JSON.parse(panelsState.value) as PanelsState;

      historyHidden.value = openedPanels.historyHidden;
      helpHidden.value = openedPanels.helpHidden;
      formHidden.value = openedPanels.formHidden;
      visibleTabLabels.value = openedPanels.visibleTabLabels;

      await Vue.nextTick();

      dockRef.value.loadLayout();
    };

    const eraseLayout = () => {
      panelsState.value = null;
    };

    const visibleTabLabels = Vue.ref([] as string[]);
    Vue.watch(tabLabels, () => {
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

    const isIncomplete = Vue.computed(() => Utils.isIncomplete(currentCall.value));

    return () => (
      <div class='w-full h-full flex' ref={root}>
        <RibbonMenu groupName='Layout'>
          <span onClick={saveLayout}>
            <IconFA name='save' style={{'padding-right': '3px'}}/>
            <span> Save </span>
          </span>
          { panelsState.value && 
          <span onClick={loadLayout}>
            <IconFA name='life-ring' style={{'padding-right': '3px'}}/>
            <span> Load </span>
          </span> }
          { panelsState.value && <span onClick={eraseLayout}>
            <IconFA name='eraser' style={{'padding-right': '3px'}}/>
            <span> Erase </span>
          </span> }
        </RibbonMenu>
        <RibbonPanel>
          <IconFA
            name='play'
            tooltip='Run step'
            onClick={run} 
          />
          { !isIncomplete.value && <ComboPopup 
            caption={ui.iconFA('arrow-to-bottom')}
            items={['Excel']}
            onSelected={({item: format}) => {
              Utils.richFunctionViewReport(
                format,
                currentCall.value.func,
                currentCall.value,
                Utils.dfToViewerMapping(currentCall.value),
              ).then((blob) => 
                DG.Utils.download(`${currentCall.value.func.nqName} - 
                  ${Utils.getStartedOrNull(currentCall.value) ?? 'Not completed'}.xlsx`, blob));
            }}
          />}
          <IconFA
            name='pen'
            tooltip={formHidden.value ? 'Open inputs': 'Close inputs'}
            onClick={() => formHidden.value = !formHidden.value} 
            style={{'background-color': !formHidden.value ? 'var(--grey-1)': null}}
          />
          <IconFA
            name='chart-pie'
            tooltip='Restore output tabs'
            onClick={() => visibleTabLabels.value = [...tabLabels.value]} 
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
        <DockManager 
          layoutStorageName={`${currentCall.value.func.nqName}_layout`}
          onPanelClosed={handlePanelClose} 
          ref={dockRef}
        >
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
              class='overflow-scroll h-full'
            />: null }
          
          { !formHidden.value ?
            <div 
              class='flex flex-col p-2 overflow-scroll h-full'
              dock-spawn-dock-type='left'
              dock-spawn-dock-ratio={0.2}
              title='Inputs'
              ref={formRef}
            >
              <InputForm 
                funcCall={currentCall.value}
                onUpdate:funcCall={(call) => emit('update:funcCall', call)}
              />
              <div class='flex sticky bottom-0 justify-end'>
                <BigButton 
                  isDisabled={!props.callState?.isRunnable || props.callState?.isRunning} 
                  onClick={run}> 
                  Run 
                </BigButton>
              </div>
            </div>: null }
          
          { 
            visibleTabLabels.value
              .map((tabLabel) => ({tabLabel, tabDfProps: categoryToDfParam.value.inputs[tabLabel] ?? 
                categoryToDfParam.value.outputs[tabLabel]}))
              .map(({tabLabel, tabDfProps}) => {
                return [
                  tabDfProps.flatMap((tabProp) => {
                    const allConfigs = Utils.getPropViewers(tabProp).config;

                    return allConfigs.map((options) => (            
                      <div 
                        class='flex flex-col pl-2 h-full w-full'
                        {...{'title': `${dfBlockTitle(tabProp)} / ${options['type']}`}}
                      >
                        <Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[tabProp.name] ?? currentCall.value.outputs[tabProp.name]}
                          class='w-full'
                        />
                      </div>));
                  }),
                  <ScalarsPanel 
                    funcCall={currentCall.value} 
                    category={tabLabel}
                    {...{title: tabLabel}}
                  />];
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
