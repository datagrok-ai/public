import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {
  Viewer, InputForm,
  BigButton, IconFA,
  RibbonPanel, DockManager, MarkDown,
  ComboPopup,
  RibbonMenu,
  ifOverlapping,
} from '@datagrok-libraries/webcomponents-vue';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {useUrlSearchParams} from '@vueuse/core';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {FittingView} from '@datagrok-libraries/compute-utils/function-views/src/fitting-view';
import {SensitivityAnalysisView} from '@datagrok-libraries/compute-utils';
import {ScalarsPanel} from './ScalarsPanel';
import {BehaviorSubject, combineLatest, EMPTY, merge, Subject} from 'rxjs';
import {useSubscription, from} from '@vueuse/rxjs'
import {catchError, switchMap, tap, map, debounceTime, withLatestFrom, share} from 'rxjs/operators';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';

type PanelsState = {
  historyHidden: boolean,
  helpHidden: boolean,
  formHidden: boolean,
  visibleTabLabels: string[],
  layout: string,
};

type ViewersData = Record<string, DG.Viewer | undefined>;

function runViewersHook(viewersHook: ViewersHook | undefined, io: string, viewers: ViewersData, meta: any) {
  if (viewersHook) {
    for (const [type, viewer] of Object.entries(viewers)) {
      if (viewer)
        viewersHook(io, type, viewer, meta);
    }
  }
}

const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';

type TabContent = Map<string,
  {type: 'dataframe', dfProp: DG.Property, config: Record<string, string | boolean> } |
  {type: 'scalars', scalarProps: DG.Property[]}
>;

const tabToProperties = (func: DG.Func) => {
  const map = {
    inputs: new Map() as TabContent,
    outputs: new Map() as TabContent,
  };

  const processDf = (dfProp: DG.Property) => {
    const dfViewers = Utils.getPropViewers(dfProp).config;
    if (dfViewers.length === 0) return;

    dfViewers.forEach((dfViewer) => {
      const dfNameWithViewer = `${dfBlockTitle(dfProp)} / ${dfViewer['type']}`;

      const tabLabel = dfProp.category === 'Misc' ?
        dfNameWithViewer: `${dfProp.category}: ${dfNameWithViewer}`;

      map.inputs.set(tabLabel, {type: 'dataframe', dfProp: dfProp, config: dfViewer});
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

      const categoryProps = map.outputs.get(category);
      if (categoryProps && categoryProps.type === 'scalars')
        categoryProps.scalarProps.push(outputProp);
      else
        map.outputs.set(category, {type: 'scalars', scalarProps: [outputProp]});
    });

  return map;
};

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
    callMeta: {
      type: Object as Vue.PropType<Record<string, BehaviorSubject<any>>>,
    },
    isTreeLocked: {
      type: Boolean,
      default: false,
    },
    isReadonly: {
      type: Boolean,
    },
    historyEnabled: {
      type: Boolean,
      default: false,
    },
    viewersHook: {
      type: Function as Vue.PropType<ViewersHook>,
    }
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
    'runClicked': () => {},
  },
  methods: {
    savePersonalState: () => {},
    loadPersonalLayout: () => {}
  },
  setup(props, {emit, expose}) {
    const currentCall = Vue.computed(() => props.funcCall);

    const tabToPropertiesMap = Vue.computed(() => tabToProperties(currentCall.value.func));

    const tabLabels = Vue.computed(() => {
      return [
        ...tabToPropertiesMap.value.inputs.keys(),
        ...tabToPropertiesMap.value.outputs.keys(),
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
    const personalPanelsStorage = (call: DG.FuncCall) => `${call.func.nqName}_personal_state`;

    const getCurrentState = () => {
      const layout = dockRef?.value?.getLayout();
      if (!layout) return null;

      return JSON.stringify({
        historyHidden: historyHidden.value,
        helpHidden: helpHidden.value,
        formHidden: formHidden.value,
        visibleTabLabels: visibleTabLabels.value,
        layout,
      });
    };

    const getSavedPersonalState = (call: DG.FuncCall): PanelsState | null => {
      const item = localStorage.getItem(personalPanelsStorage(call));
      return item ? JSON.parse(item): null;
    };

    const savePersonalState = (call: DG.FuncCall = currentCall.value) => {
      if (!dockInited.value) return;

      const state = getCurrentState();
      if (state) localStorage.setItem(personalPanelsStorage(call), state);
    };

    const removeSavedPersonalState = () => {
      localStorage.removeItem(personalPanelsStorage(currentCall.value));
    };

    let intelligentLayout = true;
    const loadPersonalLayout = async () => {
      const personalState = getSavedPersonalState(currentCall.value);
      if (!dockRef.value || !personalState || !dockInited.value) return;

      intelligentLayout = false;

      historyHidden.value = personalState.historyHidden;
      helpHidden.value = personalState.helpHidden;
      formHidden.value = personalState.formHidden;
      visibleTabLabels.value = personalState.visibleTabLabels;

      await Vue.nextTick();

      await dockRef.value.useLayout(personalState.layout);

      intelligentLayout = true;
    };

    expose({
      savePersonalState,
      loadPersonalLayout,
    })

    const dockInited = Vue.ref(false);

    const visibleTabLabels = Vue.ref([] as string[]);

    const helpText = Vue.ref(null as null | string);

    Vue.watch(currentCall, async (_, oldCall) => {
      Utils.getContextHelp(currentCall.value.func).then((loadedHelp) => {
        helpText.value = loadedHelp ?? null;
      });

      if (oldCall) savePersonalState(oldCall);
      visibleTabLabels.value = [...tabLabels.value];
    }, {immediate: true});

    Vue.watch(currentCall, async () => {
      if (dockInited.value) await loadPersonalLayout();
    }, {flush: 'post'});

    const metaStates = Vue.computed(() => props.callMeta);
    const viewersHook = Vue.computed(() => props.viewersHook);

    const metaStates$ = from(metaStates, { immediate: true });
    const viewersHook$ = from(viewersHook, { immediate: true });
    const currentCall$ = from(currentCall, { immediate: true });

    const viewerStates$ = currentCall$.pipe(
      map((call) => {
        const nextStates: Record<string, BehaviorSubject<ViewersData>> = {};
        const states = [...Object.keys(call.inputs), ...Object.keys(call.outputs)];
        for (const io of states) {
          const viewers$ = new BehaviorSubject<ViewersData>({});
          nextStates[io] = viewers$;
        }
        return nextStates;
      }),
      share(),
    );

    const viewerChanges$ = new Subject<readonly [DG.Viewer | undefined,  string, string]>();
    const viewersSub = viewerChanges$.pipe(
      withLatestFrom(viewerStates$),
      tap(([[viewer, io, type], viewerStates]) => {
        const viewers$ = viewerStates[io];
        if (viewers$) {
          const viewers = viewers$.value;
          viewers$.next({...viewers, [type]: viewer});
        }
      })
    ).subscribe();
    useSubscription(viewersSub);

    const viewHooksSub = combineLatest([viewersHook$, metaStates$, viewerStates$]).pipe(
      debounceTime(0),
      switchMap(([viewersHook, metaStates, viewerStates]) => {
        const hooksDeps: Record<string, { meta$?: BehaviorSubject<any>, viewers$?: BehaviorSubject<ViewersData> }> = {};
        for (const [io, meta$] of Object.entries(metaStates ?? {})) {
          if (!hooksDeps[io])
            hooksDeps[io] = {};
          hooksDeps[io].meta$ = meta$;
        }
        for (const [io, viewers$] of Object.entries(viewerStates ?? {})) {
          if (!hooksDeps[io])
            hooksDeps[io] = {};
          hooksDeps[io].viewers$ = viewers$;
        }
        const hookRuns = [...Object.entries(hooksDeps)].map(([io, {meta$, viewers$}]) => {
          if (!meta$ || !viewers$)
            return EMPTY;
          return combineLatest([viewers$, meta$]).pipe(
            tap(([viewers, meta]) => runViewersHook(viewersHook, io, viewers, meta)),
            catchError((e) => {
              grok.shell.error(e);
              console.error(e);
              return EMPTY;
            })
          )
        })
        return merge(...hookRuns);
      }),
    ).subscribe();
    useSubscription(viewHooksSub);

    const setViewerRef = (viewer: DG.Viewer | undefined, io: string, type: string) => {
      viewerChanges$.next([viewer, io, type] as const);
    }

    Vue.onBeforeUnmount(() => {
      savePersonalState(currentCall.value);
    });

    const isIncomplete = Vue.computed(() => props.callState?.isOutputOutdated);
    const isRunning = Vue.computed(() => props.callState?.isRunning);
    const isRunnable = Vue.computed(() => props.callState?.isRunnable);
    const currentFunc = Vue.computed(() => currentCall.value.func);

    const features = Vue.computed(() => Utils.getFeatures(currentFunc.value));
    const isSAenabled = Vue.computed(() => Utils.getFeature(features.value, 'sens-analysis', false));
    const isExportEnabled = Vue.computed(() => Utils.getFeature(features.value, 'export', true));
    const isFittingEnabled = Vue.computed(() => Utils.getFeature(features.value, 'fitting', false));

    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    return () => {
      let lastCardLabel = null as string | null;
      let scalarCardCount = 0;

      return (
        <div class='w-full h-full flex' ref={root}>
          <RibbonMenu groupName='Layout'>
            <span onClick={loadPersonalLayout}>
              <IconFA name='life-ring' style={{'padding-right': '3px'}}/>
              <span> Load personal </span>
            </span>
            <span onClick={removeSavedPersonalState}>
              <IconFA name='eraser' style={{'padding-right': '3px'}}/>
              <span> Remove personal </span>
            </span>
          </RibbonMenu>
          <RibbonMenu groupName='Panels'>
            <span
              onClick={() => formHidden.value = !formHidden.value}
              class={'flex justify-between w-full'}
            >
              <div> <IconFA name='pen' style={menuIconStyle}/> Show inputs </div>
              { !formHidden.value && <IconFA name='check'/>}
            </span>
            <span
              onClick={() => visibleTabLabels.value = [...tabLabels.value]}
              class={'flex justify-between'}
            >
              <div> <IconFA name='chart-pie'
                style={menuIconStyle}/> Show output tabs </div>
              { visibleTabLabels.value.length === tabLabels.value.length && <IconFA name='check'/>}
            </span>
            { hasContextHelp.value && <span
              onClick={() => helpHidden.value = !helpHidden.value}
              class={'flex justify-between'}
            >
              <div> <IconFA name='question' style={menuIconStyle}/> Show help </div>
              { !helpHidden.value && <IconFA name='check'/>}
            </span> }
            { props.historyEnabled && <span
              onClick={() => historyHidden.value = !historyHidden.value}
              class={'flex justify-between'}
            >
              <div> <IconFA name='history' style={menuIconStyle}/> Show history </div>
              { !historyHidden.value && <IconFA name='check'/>}
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
              name='question'
              tooltip={ helpHidden.value ? 'Open help panel' : 'Close help panel' }
              onClick={() => helpHidden.value = !helpHidden.value}
              style={{'background-color': !helpHidden.value ? 'var(--grey-1)': null}}
            /> }
            { props.historyEnabled && <IconFA
              name='history'
              tooltip='Open history panel'
              onClick={() => historyHidden.value = !historyHidden.value}
              style={{'background-color': !historyHidden.value ? 'var(--grey-1)': null}}
            /> }
          </RibbonPanel>
          <DockManager
            onPanelClosed={handlePanelClose}
            onUpdate:activePanelTitle={handleActivePanelChanged}
            onInitFinished={() => dockInited.value = true}
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

            { !formHidden.value &&
              <div
                class='flex flex-col p-2 overflow-scroll h-full'
                dock-spawn-dock-type='left'
                dock-spawn-dock-ratio={0.2}
                dock-spawn-title='Inputs'
                ref={formRef}
              >
                {
                  Vue.withDirectives(<InputForm
                    funcCall={currentCall.value}
                  />, [[ifOverlapping, isRunning.value, 'Recalculating...']])
                }
                <div class='flex sticky bottom-0 justify-end'>
                  <BigButton
                    isDisabled={!isRunnable.value || isRunning.value || props.isTreeLocked || props.isReadonly}
                    onClick={run}>
                  Run
                  </BigButton>
                </div>
              </div> }

            {
              visibleTabLabels.value
                .map((tabLabel) => ({tabLabel, tabContent: tabToPropertiesMap.value.inputs.get(tabLabel) ??
                tabToPropertiesMap.value.outputs.get(tabLabel)!}))
                .map(({tabLabel, tabContent}) => {
                  if (tabContent.type === 'dataframe') {
                    const options = tabContent.config;
                    const dfProp = tabContent.dfProp;
                    return <div
                      class='flex flex-col pl-2 h-full w-full'
                      dock-spawn-title={tabLabel}
                      key={tabLabel}
                    >
                      {
                        Vue.withDirectives(<Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[dfProp.name] ?? currentCall.value.outputs[dfProp.name]}
                          class='w-full'
                          onViewerChanged={(v) => setViewerRef(v, dfProp.name, options['type'] as string)}
                        />, [[ifOverlapping, isRunning.value, 'Recalculating...']])
                      }
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

                    return Vue.withDirectives(panel, [[ifOverlapping, isRunning.value, 'Recalculating...']]);
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
