import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BigButton, Button, InputForm} from '@datagrok-libraries/webcomponents-vue';
import {defineComponent, KeepAlive, onUnmounted, reactive, shallowRef} from 'vue';
import {Driver} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/Driver';
import {useObservable, useSubscription} from '@vueuse/rxjs';
import {isFuncCallState, isParallelPipelineState, 
  isSequentialPipelineState, PipelineState} from 
  '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {map, switchMap} from 'rxjs/operators';
import {merge} from 'rxjs';
import {FuncCallStateInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';

export const SimpleDriverApp = defineComponent({
  name: 'SimpleDriverApp',
  mounted() {
    console.log('SimpleDriverApp mounted');
  },
  unmounted() {
    console.log('SimpleDriverApp unmounted');
  },
  setup() {
    const driver = new Driver();
    const tree = shallowRef<PipelineState | undefined>(undefined);
    const state = reactive({callsInfo: {} as Record<string, FuncCallStateInfo | undefined>});

    useSubscription((driver.currentState$).subscribe((s) => tree.value = s));
    const isLockedRO = useObservable(driver.globalROLocked$);
    const isLockedMutations = useObservable(driver.treeMutationsLocked$);

    useSubscription(driver.currentCallsState$.pipe(
      switchMap((data) => {
        state.callsInfo = {};
        const entries = Object.entries(data).map(([name, state$]) => state$.pipe(map((s) => [name, s] as const)));
        return merge(...entries);
      }),
    ).subscribe(([k, val]) => {
      state.callsInfo[k] = val;
    }));

    onUnmounted(() => {
      driver.close();
      console.log('VuewDriverTestApp driver closed');
    });

    const initPipeline = (provider: string) => {
      driver.sendCommand({event: 'initPipeline', provider});
    };

    const runStep = (uuid: string) => {
      driver.sendCommand({event: 'runStep', uuid});
    };

    const addStep = (parentUuid: string, itemId: string, position: number) => {
      driver.sendCommand({event: 'addDynamicItem', parentUuid, itemId, position});
    };

    const removeStep = (uuid: string) => {
      driver.sendCommand({event: 'removeDynamicItem', uuid});
    };

    const getWidgetTree = (node: PipelineState, isParentDynamic = false) => {
      if (isFuncCallState(node)) {
        return node.funcCall ?
          <div key={node.uuid} style={{border: 'solid black 1px', margin: '5px', padding: '5px'}}>
            <div>FuncCall: {node.configId} { isParentDynamic ? <Button onClick={() => removeStep(node.uuid)}>remove</Button> : '' }</div>
            <InputForm funcCall={node.funcCall}></InputForm>
            <div>{state.callsInfo?.[node.uuid]?.isOutputOutdated ? 'OUTDATED' : 'ACTUAL'}</div>
            { Object.entries(node.funcCall.outputs).map(([name, val]) => <div>{ `${name}: ${val}`}</div>) }
            <Button onClick={() => runStep(node.uuid)}>run</Button></div> :
          <div style={{border: 'solid black 1px', margin: '5px'}}>LOADING</div>;
      } else {
        const isDynamic = (isParallelPipelineState(node) || isSequentialPipelineState(node));
        return (
          <div key={node.uuid} style={{border: 'solid black 3px', padding: '20px', margin: '5px'}}>
            <div>Pipeline: {node.configId} { isParentDynamic ? <Button onClick={() => removeStep(node.uuid)}>remove</Button> : '' }</div>
            <div style={{display: 'flex', flexDirection: 'row'}}>
              { node.steps.map((step) => getWidgetTree(step, isDynamic)) }
            </div>
            <div>
              { isDynamic ? node.stepTypes.map((step) => (<Button onClick={() => addStep(node.uuid, step.configId, node.steps.length)}>add {step.configId}</Button>)) : '' }
            </div>
          </div>
        );
      }
    };

    return () => (
      <KeepAlive>
        <div style={{width: '100%', height: '100%'}}>
          <BigButton onClick={() => initPipeline('LibTests:MockProvider3')}>Init Pipeline</BigButton>
          <div>{ isLockedRO.value ? 'TREE RO LOCKED' : 'TREE RO UNLOCKED' }</div>
          <div>{ isLockedMutations.value ? 'TREE MUTATIONS LOCKED' : 'TREE MUTATIONS UNLOCKED' }</div>
          {tree.value ? getWidgetTree(tree.value) : <div>NO TREE</div>}
        </div>
      </KeepAlive>
    );
  },
});
