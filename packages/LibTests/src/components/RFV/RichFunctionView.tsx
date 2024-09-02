import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {defineComponent, onMounted, PropType, ref, triggerRef, nextTick, computed, watch, shallowRef} from 'vue';
import {type ViewerT} from '@datagrok-libraries/webcomponents/src';
import {Viewer, InputForm, BigButton, Button, TabHeaderStripe, Tabs, IconFA, RibbonPanel, FoldableDialog, DockedPanel, SplitH, DockManager} from '@datagrok-libraries/webcomponents-vue/src';
import './RichFunctionView.css';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {History} from '../History/History';

declare global {
  namespace JSX {
    interface IntrinsicElements {
      'dg-viewer': ViewerT
    }
  }
}

export const ScalarTable = defineComponent({
  name: 'ScalarTable',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
    category: {
      type: String,
      required: true,
    },
  },
  setup(props) {
    const categoryScalars = computed(() => {
      return [
        ...props.funcCall.func.inputs
          .filter((prop) => prop.category === props.category),
        ...props.funcCall.func.outputs
          .filter((prop) => 
            prop.category === props.category ||
            (['Misc', 'Output'].includes(prop.category) && props.category === 'Output'),
          ),
      ].filter((prop) => 
        prop.propertyType !== DG.TYPE.DATA_FRAME && 
        prop.propertyType !== DG.TYPE.GRAPHICS,
      );
    });

    return () => 
      <table>
        { 
          categoryScalars.value.map((prop) => {
            const precision = prop.options.precision;

            const scalarValue = precision && 
                prop.propertyType === DG.TYPE.FLOAT && props.funcCall.outputs[prop.name] ?
              props.funcCall.outputs[prop.name].toPrecision(precision):
              props.funcCall.outputs[prop.name];
            const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;
            
            return <tr> 
              <td>
                { prop.caption ?? prop.name }{units}
              </td> 
              <td>
                { scalarValue ?? '[No value]' }
              </td> 
            </tr>;
          }) 
        }
      </table>;
  },
});

export const RichFunctionView = defineComponent({
  name: 'RichFunctionView',
  props: {
    funcCall: {
      type: Object as PropType<DG.FuncCall>,
      required: true,
    },
  },
  emits: {
    'update:funcCall': (call: DG.FuncCall) => call,
  },
  setup(props, {emit}) {
    const currentCall = computed(() => Utils.deepCopy(props.funcCall));

    const categoryToDfParam = computed(() => Utils.categoryToDfParamMap(currentCall.value.func));

    const tabLabels = computed(() => {
      return [
        ...Object.keys(categoryToDfParam.value.inputs),
        ...Object.keys(categoryToDfParam.value.outputs),
      ];
    });

    const selectedIdx = ref(0);

    const run = async () => {
      currentCall.value.newId();
      await currentCall.value.call();
      emit('update:funcCall', currentCall.value);
    };

    const formHidden = ref(false);
    const historyHidden = ref(true);

    const hasContextHelp = computed(() => Utils.hasContextHelp(currentCall.value.func));

    const dfBlockTitle = (dfProp: DG.Property) => dfProp.options['caption'] ?? dfProp.name ?? ' ';
          
    return () => (
      <div class='w-full h-full flex absolute top-0'>
        <RibbonPanel>
          { hasContextHelp.value && <IconFA 
            name='info' 
            tooltip='Open help panel' 
            onClick={async () => Utils.showHelpWithDelay((await Utils.getContextHelp(currentCall.value.func))!)}
          /> }
          <IconFA 
            name='history' 
            tooltip='Open history panel' 
            onClick={() => historyHidden.value = !historyHidden.value}
          />
        </RibbonPanel>
        <DockManager>
          { !historyHidden.value ? 
            <History 
              func={currentCall.value.func}
              showActions
              showBatchActions
              isHistory
              onRunChosen={(chosenCall) => emit('update:funcCall', chosenCall)}
            />: null }
          <div class='flex ui-div'>
            { !formHidden.value ?
              <div class='flex flex-col p-2 w-full'>
                <InputForm funcCall={currentCall.value}/>
                <div class='flex sticky bottom-0'>
                  <BigButton onClick={run}> Run </BigButton>
                </div>
              </div>: null }
            <IconFA 
              name={ formHidden.value ? 'chevron-right' : 'chevron-left' }
              tooltip={ formHidden.value ? 'Open form': 'Close form' }
              onClick={() => formHidden.value = !formHidden.value} 
              class='self-center p-2'
            />
          </div>
          <Tabs 
            items={tabLabels.value.map((label) => ({label}))} 
            selected={selectedIdx.value} 
            onUpdate:selected={(v) => selectedIdx.value = v}
            class='w-full p-1'
          >
            {{
              default: () => 
                tabLabels.value.map((tabLabel) => categoryToDfParam.value.inputs[tabLabel] ?? 
                categoryToDfParam.value.outputs[tabLabel])            
                  .flatMap((tabProps) => tabProps.map((prop) => ({prop, ...Utils.getPropViewers(prop)})))
                  .map(({prop, name, config: allConfigs}) => {
                    return [allConfigs.map((options) => 
                      <div class='flex flex-col h-1/2'>
                        <h2> { dfBlockTitle(prop) } </h2>
                        <Viewer
                          type={options['type'] as string}
                          options={options}
                          dataFrame={currentCall.value.inputs[name] ?? currentCall.value.outputs[name]}
                          class='w-full' 
                        />
                      </div>, 
                    ),
                    <ScalarTable 
                      funcCall={currentCall.value} 
                      category={tabLabels.value[selectedIdx.value]}
                    />];
                  }),
            }}
          </Tabs>
        </DockManager>
      </div>
    );
  },
});
