import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import { IconFA } from '@datagrok-libraries/webcomponents-vue';

const DEFAULT_FLOAT_PRECISION = 2;

export const ScalarsPanel = Vue.defineComponent({
  name: 'ScalarsPanel',
  props: {
    funcCall: {
      type: Object as Vue.PropType<DG.FuncCall>,
      required: true,
    },
    categoryScalars: {
      type: Array as Vue.PropType<DG.Property[]>,
      required: true,
    },
  },
  setup(props) {
    const getContent = (prop: DG.Property) => {
      const precision = prop.options.precision;

      const scalarValue = props.funcCall.outputs[prop.name];
      const formattedScalarValue = prop.propertyType === DG.TYPE.FLOAT && scalarValue ?
        precision ? scalarValue.toPrecision(precision): scalarValue.toFixed(DEFAULT_FLOAT_PRECISION):
        scalarValue;
      const units = prop.options['units'] ? ` [${prop.options['units']}]`: ``;

      return [formattedScalarValue, units];
    };

    const copyToClipboard = async (text: string) => {
      await navigator.clipboard.writeText(text);
      grok.shell.info('Value is copied to clipboard');
    };

    const hoveredIdx = Vue.ref(null as null | number);

    return () =>
      props.categoryScalars.length <= 3 ?
        <div
          class='flex flex-wrap justify-around'
        >
          { props.categoryScalars.map((prop, idx) => {
            const [scalarValue, units] = getContent(prop);

            return <div
              class='flex flex-col p-2 items-center gap-4 flex-nowrap'
              onMouseenter={() => hoveredIdx.value = idx}
              onMouseleave={() => hoveredIdx.value = null}
            >
              <div class='text-center' style={{color: 'var(--grey-4)'}}> { prop.caption ?? prop.name } </div>
              <span style={{fontSize: 'var(--font-size-large)'}}>
                { scalarValue } { units }
                <span style={{color: 'var(--grey-3)', paddingLeft: '3px'}} class='absolute'>
                  { hoveredIdx.value === idx && <IconFA
                    name='copy'
                    tooltip="Copy caption & value"
                    onClick={() => copyToClipboard(`${ prop.caption ?? prop.name } ${ scalarValue } ${ units } `)}
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
                props.categoryScalars.map((prop, idx) => {
                  const [scalarValue, units] = getContent(prop);

                  return <tr
                    onMouseenter={() => hoveredIdx.value = idx}
                    onMouseleave={() => hoveredIdx.value = null}
                  >
                    <td> <span> { prop.caption ?? prop.name } </span></td>
                    <td> <span> { units } </span></td>
                    <td> <span> { scalarValue } </span></td>
                    <td>
                      { hoveredIdx.value === idx && <IconFA
                        style={{color: 'var(--grey-3)'}}
                        name='copy'
                        tooltip="Copy caption & value"
                        onClick={() => copyToClipboard(`${ prop.caption ?? prop.name } ${ scalarValue } ${ units } `)}
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