import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useSubscription, useObservable} from '@vueuse/rxjs';
import {BehaviorSubject, merge, Observable} from 'rxjs';
import {switchMap, map} from 'rxjs/operators';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {ItemMetadata} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/view/ViewCommunication';
import {PipelineInstanceConfig, PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

function makeMergedItems<T>(input: Record<string, Observable<T>>) {
  const entries = Object.entries(input).map(([name, state$]) => state$.pipe(map((s) => [name, s] as const)));
  return merge(...entries);
}

export type DriverAPI = ReturnType<typeof useReactiveTreeDriver>;

export function useReactiveTreeDriver(
  providerFunc: Vue.Ref<string>,
  version: Vue.Ref<string | undefined>,
  instanceConfig: Vue.Ref<PipelineInstanceConfig | undefined>,
) {
  const driver = new Driver();

  const treeMutationsLocked = useObservable(driver.treeMutationsLocked$);
  const isGlobalLocked = useObservable(driver.globalROLocked$);
  const treeState = useObservable(driver.currentState$.pipe(map((x) => markFuncCallsRaw(x))));

  const currentMetaCallData = useObservable(driver.currentMetaCallData$);
  const hasNotSavedEdits = useObservable(driver.hasNotSavedEdits$);

  const logs = useObservable(driver.logger.logs$);
  const config = useObservable(driver.currentConfig$);
  const links = useObservable(driver.currentLinks$);
  const result = useObservable(driver.result$);

  const states = Vue.reactive({
    calls: {} as Record<string, FuncCallStateInfo | undefined>,
    validations: {} as Record<string, Record<string, ValidationResult> | undefined>,
    consistency: {} as Record<string, Record<string, ConsistencyInfo> | undefined>,
    meta: {} as Record<string, Record<string, BehaviorSubject<any>> | undefined>,
    descriptions: {} as Record<string, Record<string, string | string[]> | undefined>,
  });

  useSubscription(driver.nodesDescriptions$.pipe(
    switchMap((data) => {
      states.descriptions = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.descriptions[k] = val ? Vue.markRaw(val) : undefined;
  }));

  useSubscription(driver.currentCallsState$.pipe(
    switchMap((data) => {
      states.calls = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.calls[k] = val ? Vue.markRaw(val) : undefined;
  }));

  useSubscription(driver.currentValidations$.pipe(
    switchMap((data) => {
      states.validations = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.validations[k] = val ? Vue.markRaw(val) : undefined;
  }));

  useSubscription(driver.currentConsistency$.pipe(
    switchMap((data) => {
      states.consistency = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.consistency[k] = val ? Vue.markRaw(val) : undefined;
  }));

  useSubscription(driver.currentMeta$.pipe(
    switchMap((data) => {
      states.meta = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.meta[k] = val ? Vue.markRaw(val) : undefined;
  }));

  Vue.watch([() => providerFunc.value, () => version.value, () => instanceConfig.value], ([providerFunc, version, instanceConfig]) => {
    initPipeline(providerFunc, version, instanceConfig);
  });

  Vue.onMounted(() => {
    initPipeline(providerFunc.value, version.value, instanceConfig.value);
  });

  Vue.onUnmounted(() => {
    driver.close();
  });

  const initPipeline = (provider: string, version?: string, instanceConfig?: PipelineInstanceConfig) => {
    driver.sendCommand({event: 'initPipeline', provider, version, instanceConfig});
  };

  const loadPipeline = (funcCallId: string) => {
    driver.sendCommand({event: 'loadPipeline', funcCallId});
  };

  const loadAndReplaceNestedPipeline = (parentUuid: string, dbId: string, itemId: string, position: number) => {
    driver.sendCommand({event: 'loadDynamicItem', parentUuid, dbId, itemId, position, readonly: true, isReplace: true});
  };

  const savePipeline = (metaData?: ItemMetadata) => {
    driver.sendCommand({event: 'savePipeline', ...metaData});
  };

  const saveDynamicItem = (uuid:string, metaData?: ItemMetadata) => {
    driver.sendCommand({event: 'saveDynamicItem', uuid, ...metaData});
  };

  const runStep = async (uuid: string) => {
    driver.sendCommand({event: 'runStep', uuid});
  };

  const runSequence = async (startUuid: string, rerunWithConsistent?: boolean, includeNonNested?: boolean) => {
    driver.sendCommand({event: 'runSequence', startUuid, rerunWithConsistent, includeNonNested});
  };

  const runAction = (actionUuid: string) => {
    driver.sendCommand({event: 'runAction', uuid: actionUuid});
  };

  const consistencyReset = (stepUuid: string, ioName: string) => {
    driver.sendCommand({event: 'resetToConsistent', stepUuid, ioName});
  };

  const addStep = (parentUuid: string, itemId: string, position: number) => {
    driver.sendCommand({event: 'addDynamicItem', parentUuid, itemId, position});
  };

  const removeStep = (uuid: string) => {
    driver.sendCommand({event: 'removeDynamicItem', uuid});
  };

  const moveStep = (uuid: string, position: number) => {
    driver.sendCommand({event: 'moveDynamicItem', uuid, position});
  };

  const changeFuncCall = (uuid: string, call: DG.FuncCall) => {
    driver.sendCommand({event: 'updateFuncCall', stepUuid: uuid, funcCall: call});
  };

  const returnResult = () => {
    driver.sendCommand({event: 'returnResult'});
  };

  return {
    // driver,
    treeMutationsLocked,
    isGlobalLocked,
    treeState,
    currentMetaCallData,
    hasNotSavedEdits,
    states,
    logs,
    config,
    links,
    result,
    //
    loadPipeline,
    loadAndReplaceNestedPipeline,
    savePipeline,
    saveDynamicItem,
    runStep,
    runSequence,
    consistencyReset,
    runAction,
    addStep,
    removeStep,
    moveStep,
    changeFuncCall,
    returnResult,
  };
}


function markFuncCallsRaw(state?: PipelineState) {
  if (!state)
    return;
  if (state.type === 'funccall') {
    if (state.funcCall)
      Vue.markRaw(state.funcCall);
    return state;
  }
  for (const step of state.steps)
    markFuncCallsRaw(step);
  return state;
}
