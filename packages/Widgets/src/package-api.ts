import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace funcs {
  //A filter that lets you select exactly one category
  export async function radioButtonFilter(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:RadioButtonFilter', {});
  }

  //A filter that works with columns of multi-value cells (such as lists of identifiers)
  export async function multiValueFilter(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:MultiValueFilter', {});
  }

  //Shows current time
  export async function timeWidget(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:TimeWidget', {});
  }

  export async function tableSummary(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:TableSummary', {});
  }

  export async function inputDemo(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:InputDemo', {});
  }

  export async function fooInput(): Promise<any> {
    return await grok.functions.call('@datagrok/widgets:FooInput', {});
  }
}
