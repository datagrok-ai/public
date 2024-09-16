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
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {useStorage} from '@vueuse/core';

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

export const ScalarTable = Vue.defineComponent({
  name: 'ScalarTable',
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

export const RichFunctionView = Vue.defineComponent({
  name: 'RichFunctionView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
  },
  setup(props, {emit}) {
    const currentCall = Vue.computed(() => Utils.deepCopy(props.funcCall));

    const categoryToDfParam = Vue.computed(() => Utils.categoryToDfParamMap(currentCall.value.func));

    const tabLabels = Vue.computed(() => {
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

    const formHidden = Vue.ref(false);
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);

    const hasContextHelp = Vue.computed(() => Utils.hasContextHelp(currentCall.value.func));

    const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';

    const helpText = Vue.ref(null as null | string);
    Vue.watch(currentCall, async () => {
      const loadedHelp = await Utils.getContextHelp(currentCall.value.func);

      helpText.value = loadedHelp ?? null;
    }, {immediate: true});
    
    const root = Vue.ref(null as HTMLElement | null);
    const historyRef = Vue.ref(null as InstanceType<typeof History> | null);
    const helpRef = Vue.ref(null as InstanceType<typeof MarkDown> | null);
    const formRef = Vue.ref(null as HTMLElement | null);
    const dockRef = Vue.ref(null as InstanceType<typeof DockManager> | null);
    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
      if (el === formRef.value) formHidden.value = true;

      const tabIdx = visibleTabLabels.value.findIndex((label) => label === el.title);
      if (tabIdx >= 0) 
        visibleTabLabels.value.splice(tabIdx, 1);  
    };

    const save = () => {
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

    const load = async () => {
      if (!dockRef.value || !panelsState.value) return;

      const openedPanels = JSON.parse(panelsState.value) as PanelsState;

      historyHidden.value = openedPanels.historyHidden;
      helpHidden.value = openedPanels.helpHidden;
      formHidden.value = openedPanels.formHidden;
      visibleTabLabels.value = openedPanels.visibleTabLabels;

      await Vue.nextTick();

      dockRef.value.loadLayout();
    };

    const visibleTabLabels = Vue.ref([] as string[]);
    Vue.watch(tabLabels, () => {
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

    const isDataFrame = (prop: DG.Property) => (prop.propertyType === DG.TYPE.DATA_FRAME);

    const dfToViewerMapping = () => {
      const func = currentCall.value.func;

      const mapping = {} as Record<string, DG.Viewer[]>;
      Promise.all(func.inputs
        .filter((output) => isDataFrame(output))
        .map(async (p) => {
          mapping[p.name] = await Promise.all(Utils.getPropViewers(p).config
            .map((config) => configToViewer(currentCall.value.inputs[p.name], config)));

          return mapping[p.name];
        }));

      Promise.all(func.outputs
        .filter((output) => isDataFrame(output))
        .map(async (p) => {
          mapping[p.name] = await Promise.all(Utils.getPropViewers(p).config
            .map((config) => configToViewer(currentCall.value.outputs[p.name], config)));
  
          return mapping[p.name];
        }));

      return mapping;
    };

    const configToViewer = async (df: DG.DataFrame, config: Record<string, any>) => {
      const type = config['type'];
      const viewer = await df.plot.fromType(type) as DG.Viewer;
      viewer.setOptions(config);
    
      return viewer; 
    };

    return () => (
      <div class='w-full h-full flex' ref={root}>
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
          {/* <ComboPopup 
            caption={ui.iconFA('arrow-to-bottom')}
            items={['Excel']}
            onSelected={({item: format}) => {
              Utils.richFunctionViewExport(
                root.value!,
                format,
                currentCall.value.func,
                currentCall.value,
                dfToViewerMapping,
              ).then((blob) => DG.Utils.download('Test name', blob));
            }}
          /> */}
          <IconFA
            name='arrow-to-bottom'
            tooltip='Generate report'
            onClick={() => Utils.richFunctionViewExport(
              'Excel',
              currentCall.value.func,
              currentCall.value,
              dfToViewerMapping(),
            ).then((blob) => DG.Utils.download('Test name.xlsx', blob))
            }
          />
          <IconFA
            name='chart-pie'
            tooltip='Restore output tabs'
            onClick={() => visibleTabLabels.value = [...tabLabels.value]} 
          />
          <IconFA
            name='save'
            tooltip='Save the layout'
            onClick={save}
          />
          { panelsState.value && <IconFA
            name='life-ring'
            tooltip='Load the layout'
            onClick={load}
          /> }
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
              <InputForm funcCall={currentCall.value}/>
              <div class='flex sticky bottom-0'>
                <BigButton onClick={run}> Run </BigButton>
              </div>
            </div>: null }
          
          { 
            visibleTabLabels.value
              .map((tabLabel) => ({tabLabel, tabDfProps: categoryToDfParam.value.inputs[tabLabel] ?? 
                categoryToDfParam.value.outputs[tabLabel]}))            
              .map(({tabLabel, tabDfProps}) => {
                return <div 
                  class='flex flex-col overflow-scroll h-full'
                  title={tabLabel}
                  key={tabLabel}
                >
                  { tabDfProps.map((tabProp) => {
                    const allConfigs = Utils.getPropViewers(tabProp).config;

                    return allConfigs.map((options) => (            
                      <div class='flex flex-col pl-2' style={{flex: 1}}>
                        <h2> { dfBlockTitle(tabProp) } </h2>
                        <Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[tabProp.name] ?? currentCall.value.outputs[tabProp.name]}
                          class='w-full'
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
