import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {DockManager, FunctionCard, IconFA, MarkDown, RibbonMenu, tooltip} from '@datagrok-libraries/webcomponents-vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {hasAddControls, PipelineWithAdd} from '../../utils';
import {isParallelPipelineState, isSequentialPipelineState, PipelineState, ViewAction} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';


export const PipelineView = Vue.defineComponent({
  name: 'PipelineView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall | undefined>,
      required: true,
    },
    state: {
      type: Object as Vue.PropType<PipelineState>,
      required: true,
    },
    isRoot: {
      type: Boolean,
      required: true,
    },
    uuid: {
      type: String,
    },
    menuActions: {
      type:  Object as Vue.PropType<Record<string, ViewAction[]>>
    },
    buttonActions: {
      type:  Object as Vue.PropType<ViewAction[]>
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
    'proceedClicked': () => {},
    'actionRequested': (actionUuid: string) => actionUuid,
    'addNode': ({itemId, position}:{itemId: string, position: number}) => ({itemId, position})
  },
  setup(props, {emit}) {
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);
    const functionsHidden = Vue.ref(true);
    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);
    const functionsRef = Vue.shallowRef(null as HTMLElement | null);
    const state = Vue.computed(() => props.state)
    const menuActions = Vue.computed(() => props.menuActions);
    const buttonActions = Vue.computed(() => props.buttonActions);
    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    const helpText = Vue.ref(null as null | string);
    Vue.watch(() => props.funcCall, async (funcCall) => {
      if (!funcCall) {
        helpText.value = null;
        return;
      }
      const loadedHelp = await Utils.getContextHelp(funcCall.func);

      helpText.value = loadedHelp ?? null;
    }, {immediate: true});

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
    };

    const hasInnerStep = Vue.computed(() =>
      (isParallelPipelineState(state.value) || isSequentialPipelineState(state.value)) &&
      state.value.steps.length > 0
    )

    const cardsClasses = 'grok-app-card grok-gallery-grid-item-wrapper pr-4';

    return () => (
      <div class='w-full h-full flex'>
        <DockManager
          onPanelClosed={handlePanelClose}
          key={props.uuid}
        >
          { (!historyHidden.value && props.funcCall) &&
            <History
              func={props.funcCall.func}
              showActions
              showBatchActions
              isHistory
              onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
              dock-spawn-dock-type='fill'
              dock-spawn-title='History'
              dock-spawn-panel-icon='history'
              ref={historyRef}
              class='overflow-scroll h-full'
            /> }
          { (!helpHidden.value && helpText.value) &&
            <MarkDown
              markdown={helpText.value}
              dock-spawn-title='Help'
              dock-spawn-dock-type='fill'
              ref={helpRef}
            /> }
          { hasAddControls(state.value) && !functionsHidden.value && <div 
            class={'grok-gallery-grid'} 
            dock-spawn-title='Steps to add'
            dock-spawn-dock-type='fill'
            ref={functionsRef}
          >
            { state.value.stepTypes
              .map((stepType, idx) => stepType.nqName ? <FunctionCard 
                func={DG.Func.byName(stepType.nqName)}
                onClick={() => {
                  const data = state.value as PipelineWithAdd;
                  emit('addNode', {
                    itemId: data.stepTypes[idx].configId,
                    position: data.steps.length,
                  });
                }}
              /> : <div
              onClick={() => {
                const data = state.value as PipelineWithAdd;
                emit('addNode', {
                  itemId: data.stepTypes[idx].configId,
                  position: data.steps.length,
                });
              }}
              class={'grok-gallery-grid-item-wrapper'}
              style={{cursor: 'pointer'}}
            >
              <div class={'grok-gallery-grid-item grok-scripting-script d4-flex-col d4-gallery-card entity-script'}>
                <div class={'d4-flex-col'}>
                  <span class={'d4-link-label'}>
                    <label class={'grok-gallery-grid-item-title'}>
                      {stepType.friendlyName ?? stepType.configId}
                    </label>
                  </span>
                </div>
              </div>
            </div>)
            }
            </div>
          }
          { menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
            <RibbonMenu groupName={category}>
              {
                actions.map((action) => Vue.withDirectives(
                  <span onClick={() => emit('actionRequested', action.uuid)}>
                    <div> { action.icon && <IconFA name={action.icon} style={menuIconStyle}/> } { action.friendlyName ?? action.uuid } </div>
                  </span>, [[tooltip, action.description]]))
              }
            </RibbonMenu>)
          }
          <div
            dock-spawn-hide-close-button
            dock-spawn-dock-type='left'
            dock-spawn-dock-ratio={0.35}
            dock-spawn-title='Actions'
            {...{mode: 'Card'}}
            class={'p-2'}
            style={{overflow: 'auto'}}
          >
            <div
             style={{minWidth: '200px'}}
            >
              <span>
              This is a sequence of steps. You may:
              </span>

              <div class={'grok-gallery-grid'}>
                { props.funcCall &&
                  <div
                    class={cardsClasses}
                    onClick={() => historyHidden.value = false}
                  >
                    <IconFA name='history' class={'d4-picture'} />
                    <div> Load completed run </div>
                  </div> }
                { helpText.value &&
                  <div
                    class={cardsClasses}
                    onClick={() => helpHidden.value = false}
                  >
                    <IconFA name='info' class={'d4-picture'} />
                    <div> Review the docs </div>
                  </div> }
                { buttonActions.value?.map((action) => Vue.withDirectives(
                  <div
                    class={cardsClasses}
                    onClick={() => emit('actionRequested', action.uuid)}
                  >
                    <IconFA name={action.icon ?? 'circle'} class={'d4-picture'} />
                    <div> { action.friendlyName ?? action.uuid } </div>
                  </div>, [[tooltip, action.description]]))
                }
              </div>

              { hasAddControls(state.value) && <div
                class={cardsClasses}
                onClick={() => functionsHidden.value = false}
              >
                <IconFA name='plus' class={'d4-picture'} />
                <div> Choose a step to add </div>
              </div> }
              { hasInnerStep.value && <div
                class={cardsClasses}
                onClick={() => emit('proceedClicked')}
              >
                <IconFA name='plane-departure' class={'d4-picture'} />
                <div> Proceed to the sequence's first step </div>
              </div> }
            </div>
          </div>
        </DockManager>
      </div>
    );
  },
});
