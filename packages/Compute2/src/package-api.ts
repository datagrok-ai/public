import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function init(): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:Init', {});
  }

  export async function customFunctionViewEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:CustomFunctionViewEditor', { call });
  }

  export async function richFunctionViewEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:RichFunctionViewEditor', { call });
  }

  export async function treeWizardEditor(call: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TreeWizardEditor', { call });
  }

  export async function startWorkflow(nqName: string, version: string, instanceConfig: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:StartWorkflow', { nqName, version, instanceConfig });
  }

  export async function viewerTestApp(): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:ViewerTestApp', {});
  }

  export async function formTestApp(): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:FormTestApp', {});
  }

  export async function historyTestApp(): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:HistoryTestApp', {});
  }

  export async function mockPipeline1(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:MockPipeline1', { params });
  }

  export async function mockPipeline2(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:MockPipeline2', { params });
  }

  export async function testAdd2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestAdd2', { a, b });
  }

  export async function testSub2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestSub2', { a, b });
  }

  export async function testMul2(a: string, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestMul2', { a, b });
  }

  export async function testDiv2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestDiv2', { a, b });
  }

  export async function testDF1(df: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestDF1', { df });
  }

  export async function testCustomView(): Promise<any> {
    return await grok.functions.call('@datagrok/compute2:TestCustomView', {});
  }
}
