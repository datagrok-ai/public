import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {
  export async function dummy(stringInput: string, intInput: number, doubleInput: number, boolInput: boolean, choiceInput: string, tableInput: DG.DataFrame): Promise<void> {
    return await grok.functions.call('@datagrok/ui-tests:Dummy', { stringInput, intInput, doubleInput, boolInput, choiceInput, tableInput });
  }
}

export namespace funcs {
  //Viewer to test properties and others
  export async function testViewerForProperties(): Promise<any> {
    return await grok.functions.call('@datagrok/ui-tests:TestViewerForProperties', {});
  }

  //Test custom filter
  export async function testCustomFilter(): Promise<any> {
    return await grok.functions.call('@datagrok/ui-tests:TestCustomFilter', {});
  }
}
