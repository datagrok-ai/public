import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {Subject, BehaviorSubject, of, from} from 'rxjs';
import dayjs from 'dayjs';
import {RichFunctionView} from '../components/RFV/RichFunctionView';
import {getViewersHook, historyUtils, saveIsFavorite} from '@datagrok-libraries/compute-utils';
import {catchError, debounceTime, filter, switchMap, take, withLatestFrom} from 'rxjs/operators';
import {useUrlSearchParams} from '@vueuse/core';
import {EditRunMetadataDialog} from '@datagrok-libraries/compute-utils/shared-components/src/history-dialogs';
import {ViewersHook} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineConfiguration';
import {compositorOverlay} from '../directives/compositor-overlay';
import {canUseResults, pinView} from '../utils';
import {parseUrlInputs, applyUrlInputs, missingMandatoryInputs, buildInputsUrl, copyText} from '../url-inputs';
import {_package} from '../package-instance';

const RUN_DEBOUNCE_TIME = 250;
const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';

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
    view: {
      type: DG.View,
      required: true,
    },
  },
  setup(props) {
    const currentFuncCall = Vue.shallowRef(Vue.markRaw(props.funcCall));

    const viewersHook = Vue.shallowRef<ViewersHook | undefined>(undefined);

    const viewersHookMakerName = getViewersHook(currentFuncCall.value.func);

    const runRequests$ = new Subject<true>();
    const isFormValid$ = new BehaviorSubject<boolean>(false);

    const sub = runRequests$.pipe(
      debounceTime(RUN_DEBOUNCE_TIME),
      withLatestFrom(isFormValid$),
      switchMap(([, isValid]) => {
        if (!isValid || currentCallState.value.isRunning)
          return of(null);
        // Contain run() failures: an error reaching the outer subscription terminates it and kills autorun.
        return from(run()).pipe(
          catchError((err) => {
            grok.shell.error(err instanceof Error ? err.message : `${err}`);
            return of(null);
          }),
        );
      }),
    ).subscribe();

    const currentCallState = Vue.ref(
      {isRunning: false, isOutputOutdated: true, isRunnable: false, runError: undefined, pendingDependencies: []},
    );
    const overlayActive = Vue.ref(false);
    const currentView = Vue.computed(() => Vue.markRaw(props.view));

    const func = Vue.shallowRef<DG.Func | undefined>(undefined);
    const isRunningOnInput = Vue.ref<boolean>(false);

    Vue.watch(currentFuncCall, (call) => {
      func.value = call.func;
      isRunningOnInput.value = func.value.options['runOnInput'] === 'true';
    }, {immediate: true});

    const searchParams = useUrlSearchParams<{id?: string}>('history');

    const setViewName = (name: string = '') => {
      if (props.view)
        props.view.name = name;
    };

    const setViewPath = (path: string = '') => {
      if (props.view)
        props.view.path = path;
    };

    Vue.watch(searchParams, (params) => {
      setViewPath(params.id ? `?id=${params.id}` : '?');
    });

    Vue.watch(currentFuncCall, async (fc) => {
      const isOutputOutdated = JSON.parse(fc.options[OUTPUT_OUTDATED_PATH] ?? 'true');
      currentCallState.value = {...currentCallState.value, isOutputOutdated};
      searchParams.id = fc.author ? fc.id : undefined;

      const title = fc?.options?.['title'];
      const modelName = fc?.func?.friendlyName ?? fc?.func?.name;

      if (title) setViewName(`${modelName} - ${title}`);
      else setViewName(modelName);
    }, {immediate: true});

    const formReplaced$ = new BehaviorSubject<DG.InputForm | undefined>(undefined);
    // dataframe/file inputs that came from URL entity ids, kept for link export
    const urlEntityIds = new Map<string, string>();

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
      const {patch, entityIds, warnings} = await parseUrlInputs(call, urlParams);
      for (const warning of warnings)
        grok.shell.warning(warning);
      for (const [name, id] of entityIds)
        urlEntityIds.set(name, id);
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

      if (loadingId) {
        const fc = await historyUtils.loadRun(loadingId);
        currentFuncCall.value = Vue.markRaw(fc);
        return;
      }

      if ([...startUrl.searchParams.keys()].length > 0)
        await applyUrlInputsFlow(startUrl.searchParams);
    }, {immediate: true});

    const copyUrlWithInputs = async () => {
      const {url, skipped} = buildInputsUrl(currentFuncCall.value, urlEntityIds);
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
        updateCallState({isOutputOutdated: false, isRunning: false});
      } finally {
        if (currentCallState.value.isRunning)
          updateCallState({isRunning: false});
        overlayActive.value = false;
      }
    };

    const onInputChanged = () => {
      pinView(props.view);
      currentFuncCall.value.options[OUTPUT_OUTDATED_PATH] = 'true';
      updateCallState({isOutputOutdated: true});
      searchParams.id = undefined;
      if (isRunningOnInput.value)
        runRequests$.next(true);
    };

    const saveRun = async () => {
      // Invoked by RichFunctionView's shared save-to-history icon (onSaveToHistory). Block
      // saving a stale/in-flight run with a shell message.
      if (!canUseResults(currentCallState.value, 'saving'))
        return;
      const dialog = new EditRunMetadataDialog({
        title: currentFuncCall.value.options['title'] ?? '',
        description: currentFuncCall.value.options['description'] ?? '',
        tags: currentFuncCall.value.options['tags'] ?? [],
      });
      dialog.onMetadataEdit.pipe(take(1)).subscribe(async (editOptions) => {
        currentFuncCall.value.options['title'] = editOptions.title;
        currentFuncCall.value.options['description'] = editOptions.description;
        currentFuncCall.value.options['tags'] = editOptions.tags;
        currentFuncCall.value.newId();
        await historyUtils.saveRun(currentFuncCall.value);
        await saveIsFavorite(currentFuncCall.value, !!editOptions.isFavorite);
        // saveRun persists the run under currentFuncCall's id (set by newId() above), so map that
        // id straight to the URL. Set it here rather than via triggerRef -> the currentFuncCall
        // watcher, whose `fc.author` gate would clear it (a just-saved call has no author yet).
        const fc = currentFuncCall.value;
        const modelName = fc.func?.friendlyName ?? fc.func?.name;
        setViewName(fc.options['title'] ? `${modelName} - ${fc.options['title']}` : modelName);
        searchParams.id = fc.id;
      });
      dialog.show({center: true, width: 500});
    };

    const onUpdateForm = () => {
      if (isRunningOnInput.value && currentCallState.value.isOutputOutdated)
        runRequests$.next(true);
    };

    // Optional artifact publishing, double-gated: the ArtifactAlignment package must be
    // installed (its dialog function resolves) AND the enableArtifactPublishing package
    // setting must be on.
    const publishFunc = _package.settings?.['enableArtifactPublishing'] === true ?
      (DG.Func.find({name: 'publishWorkflowRunDialog'})[0] ?? null) : null;

    const publishRun = async () => {
      if (!canUseResults(currentCallState.value, 'publishing'))
        return;
      const fc = currentFuncCall.value;
      await publishFunc!.prepare({
        sourceCall: fc,
        defaultName: fc.options['title'] ?? fc.func?.friendlyName ?? fc.func?.name,
      }).call();
    };

    Vue.onUnmounted(() => {
      sub.unsubscribe();
    });

    // TODO: better async handling
    if (viewersHookMakerName) {
      const hookMaker = DG.Func.byName(viewersHookMakerName);
      hookMaker.apply().then((hook) => viewersHook.value = hook)
    }

    return () => (
      Vue.withDirectives(<div class='w-full h-full flex'>
        <RichFunctionView
          funcCall={currentFuncCall.value}
          callState={currentCallState.value}
          viewersHook={viewersHook.value}
          onUpdate:funcCall={(fc) => currentFuncCall.value = Vue.markRaw(fc)}
          onRunClicked={() => runRequests$.next(true)}
          onFormReplaced={(form) => {
            formReplaced$.next(form);
            onUpdateForm();
          }}
          urlExportHandler={copyUrlWithInputs}
          onFormValidationChanged={(val) => isFormValid$.next(val)}
          onFormInputChanged={onInputChanged}
          onSaveToHistory={() => saveRun()}
          showPublish={publishFunc != null}
          onPublishRun={() => publishRun()}
          historyEnabled={true}
          localValidation={true}
          skipInit={false}
          showRunButton={!isRunningOnInput.value}
          keepExportsVisible={isRunningOnInput.value}
          view={currentView.value}
        />
      </div>, [[compositorOverlay, overlayActive.value]])
    );
  },
});
