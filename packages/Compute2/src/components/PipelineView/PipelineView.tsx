import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {DockManager, IconFA, MarkDown, RibbonMenu, tooltip} from '@datagrok-libraries/webcomponents-vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';
import {hasAddControls, PipelineWithAdd} from '../../utils';
import {isFuncCallState, PipelineState, ViewAction} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';


export const PipelineView = Vue.defineComponent({
  name: 'PipelineView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: false,
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
      type: Object as Vue.PropType<Record<string, ViewAction[]>>,
    },
    buttonActions: {
      type: Object as Vue.PropType<ViewAction[]>,
    },
    view: {
      type: DG.ViewBase,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (_call: DG.FuncCall) => true,
    'proceedClicked': () => true,
    'actionRequested': (_actionUuid: string) => true,
    'addNode': (_data: {itemId: string, position: number}) => true,
  },
  setup(props, {emit}) {
    Vue.onRenderTriggered((event) => {
      console.log('PipelineView onRenderTriggered', event);
    });

    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);
    const functionsHidden = Vue.ref(true);

    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);
    const functionsRef = Vue.shallowRef(null as HTMLElement | null);

    const state = Vue.computed(() => props.state);
    const currentCall = Vue.computed(() => props.funcCall ? Vue.markRaw(props.funcCall) : undefined);

    const menuActions = Vue.computed(() => props.menuActions);
    const buttonActions = Vue.computed(() => props.buttonActions);
    const currentView = Vue.computed(() => Vue.markRaw(props.view));
    const isRoot = Vue.computed(() => props.isRoot);

    const helpText = Vue.ref<string | null | undefined>(null);
    const hasContextHelp = Vue.ref(false);

    const showHelp = async (func?: DG.Func) => {
      if (!func)
        helpText.value = null;
      else {
        helpText.value = await Utils.getContextHelp(func);
        helpHidden.value = false;
      }
    };

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
    };

    const hasInnerStep = Vue.ref(false);
    const name = Vue.ref('');
    const version = Vue.ref<string | undefined>(undefined);

    Vue.watch(state, (state) => {
      hasInnerStep.value = !isFuncCallState(state) && state.steps.length > 0;
      name.value = !isFuncCallState(state) ? (state.friendlyName ?? state.nqName ?? '') : '';
      version.value = !isFuncCallState(state) ? (state.version ?? '') : '';
      helpText.value = null;
    }, {immediate: true});

    Vue.watch(currentCall, (call) => {
      hasContextHelp.value = call ? Utils.hasContextHelp(call.func) : false;
    }, {immediate: true});

    const description = Vue.computed(() => {
      return `This is ${name.value} workflow${version.value ? ` version ${version.value}` : ''}. You may:`;
    });

    const cardsClasses = 'grok-app-card grok-gallery-grid-item-wrapper pr-4';
    const menuIconStyle = {width: '15px', display: 'inline-block', textAlign: 'center'};

    return () => (
      <div class='w-full h-full flex'>
        <DockManager
          onPanelClosed={handlePanelClose}
          key={props.uuid}
        >
          { (!historyHidden.value && currentCall.value) &&
            <History
              func={currentCall.value.func}
              version={!isFuncCallState(state.value) ? state.value.version : undefined}
              allowOtherVersions={isRoot.value}
              onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.7}
              dock-spawn-title='History'
              dock-spawn-panel-icon='history'
              ref={historyRef}
              class='overflow-scroll h-full'
            /> }
          { hasAddControls(state.value) && !functionsHidden.value && <div
            class={'grok-gallery-grid'}
            dock-spawn-title='Steps to add'
            dock-spawn-dock-type='right'
            dock-spawn-dock-ratio={0.7}
            ref={functionsRef}
          >
            { state.value.stepTypes
              .filter((item) => !item.disableUIAdding)
              .map((stepType, idx) => {
                const func = stepType.nqName ? Vue.markRaw(DG.Func.byName(stepType.nqName)): null;
                const language = func instanceof DG.Script ? func.language: 'javascript';
                const iconBackground = `background-image: url("/images/entities/${language}.png"); padding-right: 3px;`;

                return <div
                  class='p-2 border-solid border border-[#dbdcdf] m-4 hover:bg-[#F2F2F5]'
                  style={{width: '208px'}}
                >
                  <div class='flex flex-col'>
                    <div class='flex flex-row justify-between'>
                      <div class='flex flex-row items-center'>
                        { func ?
                          <i style={iconBackground} class={'grok-icon image-icon'}></i>:
                          <IconFA style={{'padding-right': '3px'}} name='folder-tree' />
                        }
                        <div>
                          <span class='text-[#2083d5]'> {stepType.friendlyName ?? stepType.configId} </span>
                        </div>
                      </div>
                      <div class='flex flex-row'>
                        { func?.options.help && <IconFA
                          class='d4-ribbon-item'
                          name='question-circle'
                          style={{'padding-left': '5px'}}
                          onClick={() => showHelp(func)}
                        />
                        }
                        <IconFA
                          class='d4-ribbon-item'
                          name='plus'
                          tooltip='Add step to the tree'
                          onClick={() => {
                            const data = state.value as PipelineWithAdd;
                            emit('addNode', {
                              itemId: data.stepTypes[idx].configId,
                              position: data.steps.length,
                            });
                          }}
                          style={{'padding-left': '5px'}}
                        />
                      </div>
                    </div>
                    { func && <label class='description'> {func.description} </label> }
                  </div>
                </div>;
              })
            }
          </div>
          }
          { menuActions.value && Object.entries(menuActions.value).map(([category, actions]) =>
            <RibbonMenu groupName={category} view={currentView.value}>
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
            dock-spawn-dock-type='fill'
            dock-spawn-title='Actions'
            {...{mode: 'Card'}}
            class={'p-2'}
            style={{overflow: 'auto'}}
          >
            <div
              style={{minWidth: '200px'}}
            >
              <span>
                {description.value}
              </span>

              <div class={'grok-gallery-grid'}>
                { hasContextHelp.value &&
                  <div
                    class={cardsClasses}
                    onClick={() => showHelp(currentCall.value?.func)}
                  >
                    <IconFA name='question-circle' class={'d4-picture'} />
                    <div> Show help </div>
                  </div>
                }
                { currentCall.value &&
                  <div
                    class={cardsClasses}
                    onClick={() => historyHidden.value = false}
                  >
                    <IconFA name='history' class={'d4-picture'} />
                    <div> Load completed run </div>
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
                  <div> Proceed to the first step </div>
                </div> }
              </div>
            </div>
          </div>
          { !helpHidden.value && helpText.value ?
            <div
              dock-spawn-title='Help'
              dock-spawn-dock-type='right'
              dock-spawn-dock-to='Steps to add'
              dock-spawn-dock-ratio={0.2}
              style={{overflow: 'scroll', height: '100%', paddingLeft: '5px'}}
              ref={helpRef}
            >
              <MarkDown
                markdown={helpText.value}
              />
            </div>: null
          }
        </DockManager>
      </div>
    );
  },
});
