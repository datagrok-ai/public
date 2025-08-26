import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {IconFA, ValidationIcon} from '@datagrok-libraries/webcomponents-vue';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';

export interface ScalarState {
  name: string,
  friendlyName: string,
  formattedValue: string,
  rawValue: any,
  units?: string,
}

export const ScalarsPanel = Vue.defineComponent({
  name: 'ScalarsPanel',
  props: {
    scalarsData: {
      type: Array as Vue.PropType<ScalarState[]>,
      required: true,
    },
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    }
  },
  setup(props) {
    Vue.onRenderTriggered((event) => {
      console.log('ScalarsPanel onRenderTriggered', event);
    });

    const hoveredIdx = Vue.ref(null as null | number);
    const scalarsData = Vue.computed(() => Vue.markRaw(props.scalarsData));
    const validationStates = Vue.computed(() => props.validationStates);

    const copyToClipboard = async (text: string) => {
      await navigator.clipboard.writeText(text);
      grok.shell.info('Value is copied to clipboard');
    };

    return () =>
      scalarsData.value.length <= 3 ?
        <div
          class='flex flex-wrap justify-around'
        >
          { scalarsData.value.map((prop, idx) => {
            const {formattedValue, rawValue, units, friendlyName} = prop;

            return <div
              class='flex flex-col p-2 items-center gap-4 flex-nowrap'
              onMouseenter={() => hoveredIdx.value = idx}
              onMouseleave={() => hoveredIdx.value = null}
            >
              <div class='text-center' style={{color: 'var(--grey-4)'}}> { friendlyName } </div>
              <span style={{fontSize: 'var(--font-size-large)'}}>
                { formattedValue } { units }
                <span style={{color: 'var(--grey-3)', paddingLeft: '3px'}} class='absolute'>
                  { hoveredIdx.value === idx && <IconFA
                    name='copy'
                    tooltip="Copy caption & value"
                    onClick={() => copyToClipboard(`${ friendlyName } ${ rawValue } ${ units } `)}
                  /> }
                </span>
              </span>
            </div>;
          })}
        </div> :
        <div class='h-full overflow-scroll'>
          <table class='d4-table d4-item-table d4-info-table rfv-scalar-table'>
            <tbody>
              {
                scalarsData.value.map((prop, idx) => {
                  const {formattedValue, rawValue, units, friendlyName, name} = prop;

                  return <tr
                    onMouseenter={() => hoveredIdx.value = idx}
                    onMouseleave={() => hoveredIdx.value = null}
                  >
                    <td> <span> { friendlyName } </span></td>
                    <td> <span> { units } </span></td>
                    <td> <span> { formattedValue } </span></td>
                    { validationStates.value?.[name] &&
                      <td>
                        <span>
                          <ValidationIcon validationStatus={{validation: validationStates.value?.[name]}}/>
                        </span>
                      </td>
                    }
                    <td>
                      { hoveredIdx.value === idx && <IconFA
                        style={{color: 'var(--grey-3)'}}
                        name='copy'
                        tooltip="Copy caption & value"
                        onClick={() => copyToClipboard(`${ name } ${ rawValue } ${ units } `)}
                      /> }
                    </td>
                  </tr>;
                })
              }
            </tbody>
          </table>
        </div>;
  },
});
