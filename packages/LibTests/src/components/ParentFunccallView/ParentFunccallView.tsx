import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {Button, DockManager, MarkDown} from '@datagrok-libraries/webcomponents-vue';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';

export const ParentFunccallView = Vue.defineComponent({
  name: 'ParentFunccallView',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
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
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.2}
              {...{title: 'History'}}
              ref={historyRef}
              class='overflow-scroll h-full'
            />: null }
          { !helpHidden.value && helpText.value ? 
            <MarkDown 
              markdown={helpText.value}
              {...{title: 'Help'}}
              dock-spawn-dock-type='right'
              dock-spawn-dock-ratio={0.15}
              ref={helpRef}
            /> : null 
          }
          <div>
            <span>
            This is a sequence of steps. You may:
            </span>
            <ul> 
              <li> Load a ready-to-go sequence from history 
                <Button onClick={() => historyHidden.value = false}> Choose </Button> 
              </li>
              { helpText.value && 
                <li> 
                Review the sequence documentation 
                  <Button onClick={() => helpHidden.value = false}> Open docs </Button> 
                </li> }
              <li> 
              Proceed to the sequence's first step
              
              </li>
            </ul>
          </div>
        </DockManager>
      </div>
    );
  },
});
