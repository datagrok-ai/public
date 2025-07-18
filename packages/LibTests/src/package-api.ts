import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function customInputMock(params: any): Promise<any> {
    return await grok.functions.call('LibTests:CustomInputMock', { params });
  }

  export async function rangeValidatorFactory(params: any): Promise<any> {
    return await grok.functions.call('LibTests:RangeValidatorFactory', { params });
  }

  export async function asyncValidatorDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('LibTests:AsyncValidatorDemoFactory', { params });
  }

  export async function globalValidatorDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('LibTests:GlobalValidatorDemoFactory', { params });
  }

  export async function validatorActionsDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('LibTests:ValidatorActionsDemoFactory', { params });
  }

  export async function testViewerComponent(): Promise<any> {
    return await grok.functions.call('LibTests:TestViewerComponent', {});
  }

  export async function testFromComponent(): Promise<any> {
    return await grok.functions.call('LibTests:TestFromComponent', {});
  }

  export async function testElements(): Promise<any> {
    return await grok.functions.call('LibTests:TestElements', {});
  }

  export async function testAdd2(a: number, b: number): Promise<any> {
    return await grok.functions.call('LibTests:TestAdd2', { a, b });
  }

  export async function testSub2(a: number, b: number): Promise<any> {
    return await grok.functions.call('LibTests:TestSub2', { a, b });
  }

  export async function testMul2(a: number, b: number): Promise<any> {
    return await grok.functions.call('LibTests:TestMul2', { a, b });
  }

  export async function testDiv2(a: number, b: number): Promise<any> {
    return await grok.functions.call('LibTests:TestDiv2', { a, b });
  }

  export async function testDF1(df: DG.DataFrame): Promise<any> {
    return await grok.functions.call('LibTests:TestDF1', { df });
  }

  export async function testAdd2Error(a: number, b: number): Promise<any> {
    return await grok.functions.call('LibTests:TestAdd2Error', { a, b });
  }

  export async function testMultiarg5(a: number, b: number, c: number, d: number, e: number): Promise<any> {
    return await grok.functions.call('LibTests:TestMultiarg5', { a, b, c, d, e });
  }

  export async function mockWrapper1(params: any): Promise<any> {
    return await grok.functions.call('LibTests:MockWrapper1', { params });
  }

  export async function mockWrapper2(params: any): Promise<any> {
    return await grok.functions.call('LibTests:MockWrapper2', { params });
  }

  export async function mockWrapper3(params: any): Promise<any> {
    return await grok.functions.call('LibTests:MockWrapper3', { params });
  }

  export async function mockWrapper4(params: any): Promise<any> {
    return await grok.functions.call('LibTests:MockWrapper4', { params });
  }

  export async function mockWrapper5(params: any): Promise<any> {
    return await grok.functions.call('LibTests:MockWrapper5', { params });
  }
}
