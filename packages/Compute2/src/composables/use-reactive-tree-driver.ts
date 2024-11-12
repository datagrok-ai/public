import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {useSubject, useSubscription, useObservable} from '@vueuse/rxjs';
import {BehaviorSubject, merge} from 'rxjs';
import {switchMap, map} from 'rxjs/operators';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {ConsistencyInfo, FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';

function makeMergedItems<T>(input: Record<string, BehaviorSubject<T>>) {
  const entries = Object.entries(input).map(([name, state$]) => state$.pipe(map((s) => [name, s] as const)));
  return merge(...entries);
}

export function useReactiveTreeDriver(providerFunc: Vue.Ref<string>) {
  const driver = new Driver();

  const treeMutationsLocked = useSubject(driver.treeMutationsLocked$);
  const isGlobalLocked = useSubject(driver.globalROLocked$);
  const treeState = useSubject(driver.currentState$);

  const currentMetaCallId = useSubject(driver.currentMetaCallId$);
  const hasNotSavedEdits = useSubject(driver.hasNotSavedEdits$);

  const logs = useObservable(driver.logger.logs$);
  const config = useObservable(driver.currentConfig$);
  const links = useObservable(driver.currentLinks$);

  const states = Vue.reactive({
    calls: {} as Record<string, FuncCallStateInfo | undefined>,
    validations: {} as Record<string, Record<string, ValidationResult> | undefined>,
    consistency: {} as Record<string, Record<string, ConsistencyInfo> | undefined>,
    meta: {} as Record<string, Record<string, BehaviorSubject<any>> | undefined>
  });

  useSubscription(driver.currentCallsState$.pipe(
    switchMap((data) => {
      states.calls = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.calls[k] = Object.freeze(val);
  }));

  useSubscription(driver.currentValidations$.pipe(
    switchMap((data) => {
      states.validations = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.validations[k] = Object.freeze(val);
  }));

  useSubscription(driver.currentConsistency$.pipe(
    switchMap((data) => {
      states.consistency = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.consistency[k] = Object.freeze(val);
  }));

  useSubscription(driver.currentMeta$.pipe(
    switchMap((data) => {
      states.meta = {};
      return makeMergedItems(data);
    }),
  ).subscribe(([k, val]) => {
    states.meta[k] = Object.freeze(val);
  }));

  Vue.watch(() => providerFunc.value, (providerFunc) => {
    initPipeline(providerFunc);
  });

  Vue.onMounted(() => {
    initPipeline(providerFunc.value);
  });

  Vue.onUnmounted(() => {
    driver.close();
  });

  const initPipeline = (provider: string) => {
    driver.sendCommand({event: 'initPipeline', provider});
  };

  const loadPipeline = (funcCallId: string) => {
    driver.sendCommand({event: 'loadPipeline', funcCallId});
  };

  const loadAndReplaceNestedPipeline = (parentUuid: string, dbId: string, itemId: string, position: number) => {
    driver.sendCommand({event: 'loadDynamicItem', parentUuid, dbId, itemId, position, readonly: true, isReplace: true});
  }

  const savePipeline = () => {
    driver.sendCommand({event: 'savePipeline'})
  };

  const runStep = async (uuid: string) => {
    driver.sendCommand({event: 'runStep', uuid});
  };

  const runAction = (actionUuid: string) => {
    driver.sendCommand({ event: 'runAction', uuid: actionUuid})
  }

  const addStep = (parentUuid: string, itemId: string, position: number) => {
    driver.sendCommand({event: 'addDynamicItem', parentUuid, itemId, position});
  };

  const removeStep = (uuid: string) => {
    driver.sendCommand({event: 'removeDynamicItem', uuid});
  };

  const moveStep = (uuid: string, position: number) => {
    driver.sendCommand({event: 'moveDynamicItem', uuid, position});
  };

  return {
    // driver,
    treeMutationsLocked,
    isGlobalLocked,
    treeState,
    currentMetaCallId,
    hasNotSavedEdits,
    states,
    logs,
    config,
    links,
    //
    loadPipeline,
    loadAndReplaceNestedPipeline,
    savePipeline,
    runStep,
    runAction,
    addStep,
    removeStep,
    moveStep,
  };
}
