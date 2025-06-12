import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  export async function dummy(stringInput: string, intInput: number, doubleInput: number, boolInput: boolean, choiceInput: string, tableInput: DG.DataFrame): Promise<void> {
    return await grok.functions.call('UiTests:Dummy', { stringInput, intInput, doubleInput, boolInput, choiceInput, tableInput });
  }
}

export namespace Funcs {
  //Viewer to test properties and others
  export async function testViewerForProperties(): Promise<any> {
    return await grok.functions.call('UiTests:TestViewerForProperties', {});
  }

  //Test custom filter
  export async function testCustomFilter(): Promise<any> {
    return await grok.functions.call('UiTests:TestCustomFilter', {});
  }
}
