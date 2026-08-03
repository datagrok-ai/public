import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {ValidationIcon} from '@datagrok-libraries/webcomponents-vue';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import './ScalarsPanel.css';

export interface ScalarState {
  name: string,
  friendlyName: string,
  formattedValue: string,
  rawValue: any,
  units?: string,
}

export interface ScalarsSection {
  label?: string,
  indent: number,
  scalarsData: ScalarState[],
}

export const ScalarsPanel = Vue.defineComponent({
  name: 'ScalarsPanel',
  props: {
    scalarsData: {
      type: Array as Vue.PropType<ScalarState[]>,
      required: true,
    },
    sections: {
      type: Array as Vue.PropType<ScalarsSection[]>,
      default: undefined,
    },
    validationStates: {
      type: Object as Vue.PropType<Record<string, ValidationResult>>,
    },
  },
  setup(props) {
    Vue.onRenderTriggered((event) => {
      console.log('ScalarsPanel onRenderTriggered', event);
    });

    const scalarsData = Vue.computed(() => Vue.markRaw(props.scalarsData));
    const sections = Vue.computed(() => props.sections ? Vue.markRaw(props.sections) : undefined);
    const validationStates = Vue.computed(() => props.validationStates);

    const renderScalars = (data: ScalarState[]) =>
      data.length <= 3 ?
        <div
          class='flex flex-wrap justify-around rfv2-scalar-widget'
        >
          { data.map((prop) => {
            const {formattedValue, units, friendlyName, name} = prop;

            return <div
              key={name}
              class='flex flex-col p-2 items-center gap-4 flex-nowrap'
            >
              <div class='text-center' style={{color: 'var(--grey-4)'}}> { friendlyName } </div>
              <span style={{fontSize: 'var(--font-size-large)'}}>
                { formattedValue } { units }
                <span style={{paddingLeft: '10px'}}>
                  { validationStates.value?.[name] &&
                    <td>
                      <ValidationIcon validationStatus={{validation: validationStates.value?.[name]}}/>
                    </td>
                  }
                </span>

              </span>
            </div>;
          })}
        </div> :
        <table class='d4-table d4-item-table d4-info-table rfv2-scalar-table'>
          <tbody>
            {
              data.map((prop) => {
                const {formattedValue, units, friendlyName, name} = prop;

                return <tr
                  key={name}
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
                  <td> <span> </span></td>
                </tr>;
              })
            }
          </tbody>
        </table>;

    return () => {
      if (sections.value && sections.value.length) {
        return <div class='h-full overflow-scroll p-2'>
          { sections.value.map((s, i) => (
            <div key={`${s.label ?? ''}-${s.indent}-${i}`} style={{paddingLeft: `${s.indent * 16}px`}}>
              { s.label && Vue.h(`h${Math.min(6, 2 + s.indent)}`, null, s.label) }
              { s.scalarsData.length > 0 && renderScalars(s.scalarsData) }
            </div>
          ))}
        </div>;
      }
      return scalarsData.value.length <= 3 ?
        renderScalars(scalarsData.value) :
        <div class='h-full overflow-scroll'>
          { renderScalars(scalarsData.value) }
        </div>;
    };
  },
});
