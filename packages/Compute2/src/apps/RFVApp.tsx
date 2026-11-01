import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {EMPTY, Subject, BehaviorSubject, of, from} from 'rxjs';
import dayjs from 'dayjs';
import {RichFunctionView} from '../components/RFV/RichFunctionView';
import {getRunName, getViewersHook} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {catchError, debounceTime, filter, switchMap, take, withLatestFrom} from 'rxjs/operators';
import {useUrlSearchParams} from '@vueuse/core';
import {useSubscription} from '@vueuse/rxjs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import type {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {compositorOverlay} from '../directives/compositor-overlay';
import {canUseResults, pinView} from '../utils';
import {parseUrlInputs, applyUrlInputs, missingMandatoryInputs, buildInputsUrl, copyText} from '../url-inputs';
import {defaultRunName, modelLabel} from '../components/RunStorage/run-storage';
import {openSaveRunDialog} from '../components/RunStorage/SaveRunDialog';
import {RibbonMenu, useDgView} from '@datagrok-libraries/webcomponents-vue';
import {useFeedbackItems} from '../composables/use-feedback';
import {
  OUTPUT_OUTDATED_PATH, readOutputOutdated,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {entryFromLiveCall} from '../components/RunComparison/entry-extraction';
import {
  addToComparison, snapshotName,
} from '../components/RunComparison/comparison-views';

const RUN_DEBOUNCE_TIME = 250;
const FUNCTION_HISTORY = {mode: 'function'} as const;

export const RFVApp = Vue.defineComponent({
  name: 'RFVApp',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    initialRunId: {
      type: String,
      required: false,
    },
  },
  setup(props) {
    const currentFuncCall = Vue.shallowRef(Vue.markRaw(props.funcCall));

    const viewersHook = Vue.shallowRef<ViewersHook | undefined>(undefined);

    const viewersHookMakerName = getViewersHook(currentFuncCall.value.func);

    const runRequests$ = new Subject<true>();
    const isFormValid$ = new BehaviorSubject<boolean>(false);

    useSubscription(runRequests$.pipe(
      debounceTime(RUN_DEBOUNCE_TIME),
      withLatestFrom(isFormValid$),
      switchMap(([, isValid]) => {
        if (!isValid || currentCallState.value.isRunning)
          return of(null);
        // Contain run() failures: an error reaching the outer subscription terminates it and kills autorun.
        return from(run()).pipe(
          catchError((err) => {
            updateCallState({runError: String(err)});
            grok.shell.error(err instanceof Error ? err.message : `${err}`);
            return of(null);
          }),
        );
      }),
    ).subscribe());

    const currentCallState = Vue.ref<FuncCallStateInfo>(
      {isRunning: false, isOutputOutdated: true, isRunnable: false, runError: undefined, pendingDependencies: []},
    );
    const overlayActive = Vue.ref(false);
    const currentView = useDgView();

    const func = Vue.shallowRef<DG.Func | undefined>(undefined);
    const isRunningOnInput = Vue.ref<boolean>(false);

    Vue.watch(currentFuncCall, (call) => {
      func.value = call.func;
      isRunningOnInput.value = func.value.options['runOnInput'] === 'true';
    }, {immediate: true});

    const searchParams = useUrlSearchParams<{id?: string}>('history');

    // the name of the loaded or last saved run, and whether its inputs changed since
    const savedName = Vue.ref('');
    const isEdited = Vue.ref(false);

    const setViewName = (name: string = '') => {
      currentView.name = name;
    };

    const setViewPath = (path: string = '') => {
      currentView.path = path;
    };

    Vue.watch(searchParams, (params) => {
      setViewPath(params.id ? `?id=${params.id}` : '?');
    });

    Vue.watch(currentFuncCall, async (fc) => {
      // a call that never ran has no flag yet
      const isOutputOutdated = fc.options[OUTPUT_OUTDATED_PATH] == null || readOutputOutdated(fc);
      currentCallState.value = {...currentCallState.value, isOutputOutdated};
      searchParams.id = fc.author ? fc.id : undefined;
      savedName.value = fc.author ? getRunName(fc) : '';
      isEdited.value = false;
    }, {immediate: true});

    Vue.watch([savedName, isEdited], ([name, edited]) => {
      const modelName = currentFuncCall.value.func?.friendlyName ?? currentFuncCall.value.func?.name;
      setViewName(name ? `${modelName} - ${name}${edited ? ' (edited)' : ''}` : modelName);
    }, {immediate: true});

    const formReplaced$ = new BehaviorSubject<DG.InputForm | undefined>(undefined);

    const clearUrlInputs = () => {
      for (const key of Object.keys(searchParams)) {
        if (key !== 'id')
          (searchParams as any)[key] = undefined;
      }
      setViewPath(searchParams.id ? `?id=${searchParams.id}` : '?');
    };

    // URL inputs act as overrides: entity values are pre-loaded first, then the whole
    // patch is applied in one sync block after the form is built (defaults already in)
    const applyUrlInputsFlow = async (urlParams: URLSearchParams) => {
      const call = currentFuncCall.value;
      const {patch, warnings} = await parseUrlInputs(call, urlParams);
      for (const warning of warnings)
        grok.shell.warning(warning);
      if (patch.size > 0) {
        await formReplaced$.pipe(filter((form) => form != null), take(1)).toPromise();
        if (currentFuncCall.value !== call)
          return;
        applyUrlInputs(call, patch);
        if (missingMandatoryInputs(call).length === 0)
          runRequests$.next(true);
      }
      clearUrlInputs();
    };

    // Programmatic open (OpenWorkflowRun) passes the id via call.aux, deep links via the start URL
    let initialRunId = props.initialRunId;

    Vue.watch(currentFuncCall, async () => {
      if (!initialRunId && globalThis.initialURLHandled)
        return;

      const startUrl = new URL(grok.shell.startUri);
      const loadingId = initialRunId ?? startUrl.searchParams.get('id');
      initialRunId = undefined;
      globalThis.initialURLHandled = true;

      try {
        if (loadingId) {
          const fc = await historyUtils.loadRun(loadingId);
          currentFuncCall.value = Vue.markRaw(fc);
          return;
        }

        if ([...startUrl.searchParams.keys()].length > 0)
          await applyUrlInputsFlow(startUrl.searchParams);
      } catch (e: any) {
        grok.shell.error(e);
      }
    }, {immediate: true});

    const copyUrlWithInputs = async () => {
      const {url, skipped} = buildInputsUrl(currentFuncCall.value);
      if (!await copyText(url)) {
        grok.shell.warning('Could not access the clipboard');
        return;
      }
      grok.shell.info(skipped.length > 0 ?
        `Link copied (inputs not included: ${skipped.join(', ')})` :
        'Link copied to clipboard');
    };


    // Replace the object on every transition so RichFunctionView's
    // by-reference watch on props.callState fires (it can't deep-watch).
    const updateCallState = (patch: Partial<typeof currentCallState.value>) => {
      currentCallState.value = {...currentCallState.value, ...patch};
    };

    const run = async () => {
      if (!currentFuncCall.value) return;
      overlayActive.value = true;
      updateCallState({isRunning: true});
      try {
        await currentFuncCall.value.call(undefined, undefined, {processed: true, report: false});
        currentFuncCall.value.started = dayjs();

        currentFuncCall.value.options[OUTPUT_OUTDATED_PATH] = 'false';
        updateCallState({isOutputOutdated: false, isRunning: false, runError: undefined});
      } finally {
        if (currentCallState.value.isRunning)
          updateCallState({isRunning: false});
        overlayActive.value = false;
      }
    };

    const currentSnapshot = () => {
      const fc = currentFuncCall.value;
      const {name, description, tags} = historyUtils.getRunMetadata(fc);
      const metadata = {title: name, description, tags};
      return entryFromLiveCall(fc, currentCallState.value, snapshotName(name || modelLabel(fc.func)), metadata);
    };

    const onInputChanged = () => {
      pinView(currentView);
      currentFuncCall.value.options[OUTPUT_OUTDATED_PATH] = 'true';
      updateCallState({isOutputOutdated: true});
      searchParams.id = undefined;
      isEdited.value = true;
      if (isRunningOnInput.value)
        runRequests$.next(true);
    };

    const openRequests$ = new Subject<DG.FuncCall>();
    useSubscription(openRequests$.pipe(
      // the run browser hands over the run as listed; the form needs its inputs and outputs
      switchMap((listed) => from(historyUtils.loadRun(listed.id)).pipe(
        catchError((e) => {
          grok.shell.error(`Couldn't open the run: ${e?.message ?? e}`);
          return EMPTY;
        }),
      )),
    ).subscribe((loaded) => currentFuncCall.value = Vue.markRaw(loaded)));

    const saveRunWithDialog = async (): Promise<string | null> => {
      // Invoked by RichFunctionView's save icon (onSaveToHistory). Block saving a stale/in-flight
      // run with a shell message.
      if (!canUseResults(currentCallState.value, 'saving'))
        return null;
      const fc = currentFuncCall.value;
      const {description, tags} = historyUtils.getRunMetadata(fc);
      const result = await openSaveRunDialog({
        func: fc.func, name: defaultRunName(fc.func, savedName.value), description, tags,
      });
      if (!result)
        return null;
      try {
        await historyUtils.saveRun(fc, {newId: true, metadata: result.metadata, target: result.target});
      } catch (e: any) {
        grok.shell.error(e);
        return null;
      }
      // saveRun persists the run under currentFuncCall's id (set by newId() above), so map that
      // id straight to the URL. Set it here rather than via triggerRef -> the currentFuncCall
      // watcher, whose `fc.author` gate would clear it (a just-saved call has no author yet).
      savedName.value = result.metadata.name;
      isEdited.value = false;
      searchParams.id = fc.id;
      return fc.id;
    };

    const onUpdateForm = () => {
      if (isRunningOnInput.value && currentCallState.value.isOutputOutdated)
        runRequests$.next(true);
    };

    const feedbackItems = useFeedbackItems(Vue.computed(() => currentFuncCall.value.func));

    const loadViewersHook = async () => {
      if (!viewersHookMakerName)
        return;
      try {
        viewersHook.value = await DG.Func.byName(viewersHookMakerName).apply();
      } catch (e: any) {
        grok.shell.error(e);
      }
    };
    loadViewersHook();

    return () => (
      Vue.withDirectives(<div class='w-full h-full flex'>
        <RichFunctionView
          funcCall={currentFuncCall.value}
          callState={currentCallState.value}
          viewersHook={viewersHook.value}
          onUpdate:funcCall={(fc) => currentFuncCall.value = Vue.markRaw(fc)}
          onStepRunChosen={(run) => openRequests$.next(run)}
          onRunClicked={() => runRequests$.next(true)}
          onFormReplaced={(form) => {
            formReplaced$.next(form);
            onUpdateForm();
          }}
          urlExportHandler={copyUrlWithInputs}
          onFormValidationChanged={(val) => isFormValid$.next(val)}
          onFormInputChanged={onInputChanged}
          onSaveToHistory={() => saveRunWithDialog()}
          onCompareRequested={() => addToComparison(currentSnapshot(), currentFuncCall.value.func)}
          history={FUNCTION_HISTORY}
          standalone={true}
        />
        {feedbackItems.value.length > 0 &&
          <RibbonMenu groupName='Feedback' items={feedbackItems.value}/>}
      </div>, [[compositorOverlay, overlayActive.value]])
    );
  },
});
