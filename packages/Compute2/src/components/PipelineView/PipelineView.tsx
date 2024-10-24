import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {Button, DockManager, IconFA, MarkDown} from '@datagrok-libraries/webcomponents-vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';

export const PipelineView = Vue.defineComponent({
  name: 'PipelineView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
    'proceedClicked': () => {},
  },
  setup(props, {emit}) {
    const historyHidden = Vue.ref(true);
    const helpHidden = Vue.ref(true);
    const historyRef = Vue.shallowRef(null as InstanceType<typeof History> | null);
    const helpRef = Vue.shallowRef(null as InstanceType<typeof MarkDown> | null);

    const helpText = Vue.ref(null as null | string);
    Vue.watch(() => props.funcCall, async () => {
      const loadedHelp = await Utils.getContextHelp(props.funcCall.func);

      helpText.value = loadedHelp ?? null;
    }, {immediate: true});

    const handlePanelClose = async (el: HTMLElement) => {
      if (el === historyRef.value?.$el) historyHidden.value = true;
      if (el === helpRef.value?.$el) helpHidden.value = true;
    };

    return () => (
      <div class='w-full h-full flex'>
        <DockManager 
          onPanelClosed={handlePanelClose} 
        >
          { !historyHidden.value ? 
            <History 
              func={props.funcCall.func}
              showActions
              showBatchActions
              isHistory
              onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
              dock-spawn-dock-type='fill'
              dock-spawn-title='History'
              ref={historyRef}
              class='overflow-scroll h-full'
            />: null }
          { !helpHidden.value && helpText.value ? 
            <MarkDown 
              markdown={helpText.value}
              dock-spawn-title='Help'
              dock-spawn-dock-type='fill'
              ref={helpRef}
            /> : null 
          }
          <div 
            dock-spawn-hide-close-button
            dock-spawn-dock-type='left'
            dock-spawn-dock-ratio={0.35}
            dock-spawn-title='Actions'
            {...{mode: 'Card'}}
            class={'p-2'}
            style={{minWidth: '200px'}}
          >
            <span>
            This is a sequence of steps. You may:
            </span>
             
            <div class={'grok-gallery-grid'}>
              <div 
                class='grok-app-card grok-gallery-grid-item-wrapper pr-4'
                onClick={() => historyHidden.value = false}
              > 
                <IconFA name='history' class={'d4-picture'} />
                <div> Load completed run </div>
              </div>
              { helpText.value && 
                <div 
                  class='grok-app-card grok-gallery-grid-item-wrapper pr-4'
                  onClick={() => helpHidden.value = false}
                > 
                  <IconFA name='info' class={'d4-picture'} />
                  <div> Review the docs </div>
                </div> }
              <div 
                class='grok-app-card grok-gallery-grid-item-wrapper pr-4'
                onClick={() => emit('proceedClicked')}
              > 
                <IconFA name='plane-departure' class={'d4-picture'} />
                <div> Proceed to the sequence's first step </div>
              </div>
            </div>
          </div>
        </DockManager>
      </div>
    );
  },
});
