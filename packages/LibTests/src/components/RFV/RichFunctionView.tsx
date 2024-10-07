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
import {useElementHover, useStorage, useUrlSearchParams} from '@vueuse/core';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';

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

const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';

type TabContent = Record<string, 
  {type: 'dataframe', dfProp: DG.Property, config: Record<string, string | boolean> } | 
  {type: 'scalars', scalarProps: DG.Property[]}
>;

const tabToProperties = (func: DG.Func) => {
  const map = {
    inputs: {} as TabContent,
    outputs: {} as TabContent,
  };

  const processDf = (dfProp: DG.Property) => {
    const dfViewers = Utils.getPropViewers(dfProp).config;
    if (dfViewers.length === 0) return;

    dfViewers.forEach((dfViewer) => {
      const dfNameWithViewer = `${dfBlockTitle(dfProp)} / ${dfViewer['type']}`;

      const tabLabel = dfProp.category === 'Misc' ? 
        dfNameWithViewer: `${dfProp.category}: ${dfNameWithViewer}`;

      map.inputs[tabLabel] = {type: 'dataframe', dfProp: dfProp, config: dfViewer};
    });
    return;
  };

  func.inputs
    .forEach((inputProp) => {
      if (inputProp.propertyType === DG.TYPE.DATA_FRAME) processDf(inputProp);
    });

  func.outputs
    .forEach((outputProp) => {
      if (outputProp.propertyType === DG.TYPE.DATA_FRAME) {
        processDf(outputProp);
        return;
      }
      
      const category = outputProp.category === 'Misc' ? 'Output': outputProp.category;

      if (map.outputs[category] && map.outputs[category].type === 'scalars')
        map.outputs[category].scalarProps.push(outputProp);
      else
        map.outputs[category] = {type: 'scalars', scalarProps: [outputProp]};
    });

  return map;
};

export const ScalarsPanel = Vue.defineComponent({
  name: 'ScalarsPanel',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    categoryScalars: {
      type: Array as Vue.PropType<DG.Property[]>,
      required: true,
    },
  },
  setup(props) {
    const getContent = (prop: DG.Property) => {
      const precision = prop.options.precision;

      const scalarValue = precision && 
                prop.propertyType === DG.TYPE.FLOAT && props.funcCall.outputs[prop.name] ?
        props.funcCall.outputs[prop.name].toPrecision(precision):
        props.funcCall.outputs[prop.name];
      const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

      return [scalarValue, units];
    };

    const copyToClipboard = async (text: string) => {
      await navigator.clipboard.writeText(text);
      grok.shell.info('Value is copied to clipboard');
    };

    const hoveredIdx = Vue.ref(null as null | number);
    
    return () => 
      props.categoryScalars.length <= 3 ?
        <div 
          class='flex flex-wrap justify-around' 
        >
          { props.categoryScalars.map((prop, idx) => {
            const [scalarValue, units] = getContent(prop);

            return <div 
              class='flex flex-col p-2 items-center gap-4 flex-nowrap'
              onMouseenter={() => hoveredIdx.value = idx}
              onMouseleave={() => hoveredIdx.value = null}
            > 
              <div class='text-center' style={{color: 'var(--grey-4)'}}> { prop.caption ?? prop.name } </div>
              <span style={{fontSize: 'var(--font-size-large)'}}> 
                { scalarValue } { units } 
                <span style={{color: 'var(--grey-3)', paddingLeft: '3px'}} class='absolute'> 
                  { hoveredIdx.value === idx && <IconFA 
                    name='copy' 
                    tooltip="Copy caption & value" 
                    onClick={() => copyToClipboard(`${ prop.caption ?? prop.name } ${ scalarValue } ${ units } `)}
                  /> }
                </span>
              </span>
            </div>;
          })}
        </div> :
        <div class='h-full overflow-scroll'>
          <table class='d4-table d4-item-table d4-info-table rfv-scalar-table'> 
            <tbody>
              { 
                props.categoryScalars.map((prop, idx) => { 
                  const [scalarValue, units] = getContent(prop);

                  return <tr
                    onMouseenter={() => hoveredIdx.value = idx}
                    onMouseleave={() => hoveredIdx.value = null}
                  >
                    <td> <span> { prop.caption ?? prop.name } </span></td>
                    <td> <span> { units } </span></td>
                    <td> <span> { scalarValue } </span></td>
                    <td> 
                      { hoveredIdx.value === idx && <IconFA 
                        style={{color: 'var(--grey-3)'}}
                        name='copy' 
                        tooltip="Copy caption & value" 
                        onClick={() => copyToClipboard(`${ prop.caption ?? prop.name } ${ scalarValue } ${ units } `)}
                      /> }
                    </td>
                  </tr>;
                })
              }
            </tbody>
          </table>
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

    const tabToPropertiesMap = Vue.computed(() => tabToProperties(currentCall.value.func));

    const tabLabels = Vue.computed(() => {
      return [
        ...Object.keys(tabToPropertiesMap.value.inputs),
        ...Object.keys(tabToPropertiesMap.value.outputs),
      ];
    });

    const run = async () => {
      emit('runClicked');
    };

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);

    const hasContextHelp = Vue.computed(() => Utils.hasContextHelp(currentCall.value.func));
    
    const root = Vue.ref(null as HTMLElement | null);
    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);
    const formRef = Vue.shallowRef(null as HTMLElement | null);
    const dockRef = Vue.shallowRef(null as InstanceType<typeof DockManager> | null);
    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const tabIdx = visibleTabLabels.value.findIndex((label) => label === el.getAttribute('dock-spawn-title'));
      if (tabIdx >= 0) 
        visibleTabLabels.value.splice(tabIdx, 1);  
    };

    const hashParams = useUrlSearchParams('hash-params');
    const handleActivePanelChanged = (panelTitle: string | null) => {
      hashParams.activePanel = panelTitle ?? [];
    };

    const panelsStorageName = Vue.computed(() => `${currentCall.value.func.nqName}_panels`);

    const panelsState = useStorage(
      panelsStorageName.value,
      null as null | string,
    );

    let intelligentLayout = true;
    const loadLayout = async () => {
      if (!dockRef.value || !panelsState.value) return;

      intelligentLayout = false;
      const openedPanels = JSON.parse(panelsState.value) as PanelsState;

      historyHidden.value = openedPanels.historyHidden;
      helpHidden.value = openedPanels.helpHidden;
      formHidden.value = openedPanels.formHidden;
      visibleTabLabels.value = openedPanels.visibleTabLabels;

      await Vue.nextTick();

      await dockRef.value.loadLayout();

      intelligentLayout = true;
    };

    const eraseLayout = () => {
      panelsState.value = null;
    };

    const visibleTabLabels = Vue.ref([] as string[]);
    Vue.watch(tabLabels, () => {
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

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

    const helpText = Vue.ref(null as null | string);
    Vue.watch(currentCall, async () => {
      Utils.getContextHelp(currentCall.value.func).then((loadedHelp) => {
        helpText.value = loadedHelp ?? null;
      });
      
      saveLayout();
    }, {immediate: true});

    Vue.watch(currentCall, async () => {
      await loadLayout();
    }, {'flush': 'post'});

    Vue.onBeforeUnmount(() => {
      console.log('Unmounting RFV');
    });

    const isIncomplete = Vue.computed(() => Utils.isIncomplete(currentCall.value));

    const currentFunc = Vue.computed(() => currentCall.value.func);
    const features = Vue.computed(() => Utils.getFeatures(currentFunc.value));
    const isSAenabled = Vue.computed(() => Utils.getFeature(features.value, 'sens-analysis', false));
    const isExportEnabled = Vue.computed(() => Utils.getFeature(features.value, 'export', true));
    const isFittingEnabled = Vue.computed(() => Utils.getFeature(features.value, 'fitting', false));

    return () => {
      let lastCardLabel = null as string | null;
      let scalarCardCount = 0;

      return (
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
            { panelsState.value && 
              <span onClick={eraseLayout}>
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
            { isExportEnabled.value && !isIncomplete.value && <ComboPopup 
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
            { isSAenabled.value && <IconFA 
              name='analytics'
              onClick={() => SensitivityAnalysisView.fromEmpty(currentFunc.value)}
              tooltip='Run sensitivity analysis'
            />}
            { isFittingEnabled.value && <IconFA 
              name='chart-line'
              onClick={() => FittingView.fromEmpty(currentFunc.value)}
              tooltip='Fit inputs'
            />}
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
            onUpdate:activePanelTitle={handleActivePanelChanged}
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
                dock-spawn-title='History'
                ref={historyRef}
                class='overflow-scroll h-full'
              />: null }
          
            { !formHidden.value ?
              <div 
                class='flex flex-col p-2 overflow-scroll h-full'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-title='Inputs'
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
                .map((tabLabel) => ({tabLabel, tabContent: tabToPropertiesMap.value.inputs[tabLabel] ?? 
                tabToPropertiesMap.value.outputs[tabLabel]}))
                .map(({tabLabel, tabContent}) => {
                  if (tabContent.type === 'dataframe') {
                    const options = tabContent.config;
                    const dfProp = tabContent.dfProp;
                    return <div 
                      class='flex flex-col pl-2 h-full w-full'
                      dock-spawn-title={tabLabel}
                      key={tabLabel}
                    >
                      <Viewer
                        type={options['type'] as string}
                        options={options}
                        dataFrame={currentCall.value.inputs[dfProp.name] ?? currentCall.value.outputs[dfProp.name]}
                        class='w-full'
                      />
                    </div>;
                  }

                  if (tabContent.type === 'scalars') {
                    const categoryProps = tabContent.scalarProps;
                    
                    const panel = <ScalarsPanel
                      class='h-full overflow-scroll'
                      categoryScalars={categoryProps}
                      funcCall={currentCall.value}
                      key={tabLabel}
                      dock-spawn-title={tabLabel}
                      dock-spawn-dock-to={intelligentLayout && 
                        lastCardLabel && scalarCardCount < 3 ? lastCardLabel: null
                      }
                      dock-spawn-dock-type={intelligentLayout ? (lastCardLabel ? 
                        (categoryProps.length > 3 ? 
                          'fill': (scalarCardCount < 3 ? 'right': 'down')
                        ): 'down') : null
                      }
                      dock-spawn-dock-ratio={intelligentLayout ? 
                        (lastCardLabel && scalarCardCount < 3 ? 0.5: 0.15) : null
                      }
                    />;

                    if (categoryProps.length < 3) {
                      scalarCardCount += categoryProps.length;     
                      lastCardLabel = tabLabel;
                    }
                    if (scalarCardCount >= 2) {
                      scalarCardCount = 0;
                      lastCardLabel = null;
                    }

                    return panel;
                  }
                })
            }
            { !helpHidden.value && helpText.value ? 
              <MarkDown 
                markdown={helpText.value}
                dock-spawn-title='Help'
                dock-spawn-dock-type='right'
                dock-spawn-dock-ratio={0.15}
                ref={helpRef}
              /> : null 
            }
          </DockManager>
        </div>
      );
    };
  },
});
