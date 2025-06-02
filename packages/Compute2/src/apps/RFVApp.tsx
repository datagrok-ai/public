import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {Subject, BehaviorSubject, of, from} from 'rxjs';

import {RichFunctionView} from '../components/RFV/RichFunctionView';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {debounceTime, switchMap, take, withLatestFrom} from 'rxjs/operators';
import {IconFA, RibbonPanel} from '@datagrok-libraries/webcomponents-vue';
import {useUrlSearchParams} from '@vueuse/core';
import {saveIsFavorite} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {EditRunMetadataDialog} from '@datagrok-libraries/compute-utils/shared-components/src/history-dialogs';

const RUN_DEBOUNCE_TIME = 250;
const OUTPUT_OUTDATED_PATH = 'OUTPUT_OUTDATED';

export const RFVApp = Vue.defineComponent({
  name: 'RFVApp',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    view: {
      type: DG.ViewBase,
      required: true,
    },
  },
  setup(props) {
    const runRequests$ = new Subject<true>();
    const isFormValid$ = new BehaviorSubject<boolean>(false);

    const sub = runRequests$.pipe(
      debounceTime(RUN_DEBOUNCE_TIME),
      withLatestFrom(isFormValid$),
      switchMap(([, isValid]) => {
        if (!isValid || currentCallState.value.isRunning)
          return of(null);
        return from(run());
      }),
    ).subscribe();

    const currentFuncCall = Vue.shallowRef(Vue.markRaw(props.funcCall));
    const currentCallState = Vue.ref(
      {isRunning: false, isOutputOutdated: true, isRunnable: false, runError: undefined, pendingDependencies: []},
    );
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
      currentCallState.value.isOutputOutdated = isOutputOutdated;
      searchParams.id = fc.author ? fc.id : undefined;

      const title = fc?.options?.['title'];
      const modelName = fc?.func?.friendlyName ?? fc?.func?.name;

      if (title) setViewName(`${modelName} - ${title}`);
      else setViewName(modelName);
    }, {immediate: true});

    Vue.watch(currentFuncCall, async () => {
      if (globalThis.initialURLHandled)
        return;

      globalThis.initialURLHandled = true;

      const startUrl = new URL(grok.shell.startUri);
      const loadingId = startUrl.searchParams.get('id');

      if (!loadingId)
        return;

      const fc = await historyUtils.loadRun(loadingId);
      currentFuncCall.value = Vue.markRaw(fc);
    }, {immediate: true});


    const run = async () => {
      if (!currentFuncCall.value) return;
      currentCallState.value.isRunning = true;
      try {
        await currentFuncCall.value.call();
        currentCallState.value.isOutputOutdated = false;
        currentFuncCall.value.options[OUTPUT_OUTDATED_PATH] = 'false';
      } finally {
        currentCallState.value.isRunning = false;
      }
    };

    const onInputChanged = () => {
      currentCallState.value.isOutputOutdated = true;
      currentFuncCall.value.options[OUTPUT_OUTDATED_PATH] = 'true';
      searchParams.id = undefined;
      if (isRunningOnInput.value)
        runRequests$.next(true);
    };

    const saveRun = async () => {
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
        await historyUtils.saveRun(currentFuncCall.value, false);
        await saveIsFavorite(currentFuncCall.value, !!editOptions.isFavorite);
        Vue.triggerRef(currentFuncCall);
      });
      dialog.show({center: true, width: 500});
    };

    const onUpdateForm = () => {
      if (isRunningOnInput.value && currentCallState.value.isOutputOutdated)
        runRequests$.next(true);
    };

    Vue.onUnmounted(() => {
      sub.unsubscribe();
    });

    return () => (
      <div class='w-full h-full flex'>
        <RibbonPanel view={currentView.value}>
          {!currentCallState.value.isOutputOutdated &&
            <IconFA
              name='save'
              tooltip={'Save'}
              onClick={saveRun}
            />}
        </RibbonPanel>
        <RichFunctionView
          funcCall={currentFuncCall.value}
          callState={currentCallState.value}
          onUpdate:funcCall={(fc) => currentFuncCall.value = Vue.markRaw(fc)}
          onRunClicked={() => runRequests$.next(true)}
          onFormReplaced={onUpdateForm}
          onFormValidationChanged={(val) => isFormValid$.next(val)}
          onFormInputChanged={onInputChanged}
          historyEnabled={true}
          localValidation={true}
          skipInit={false}
          showRunButton={!isRunningOnInput.value}
          view={currentView.value}
        />
      </div>
    );
  },
});
