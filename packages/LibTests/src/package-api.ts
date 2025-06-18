import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace funcs {
  export async function customInputMock(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:CustomInputMock', { params });
  }

  export async function rangeValidatorFactory(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:RangeValidatorFactory', { params });
  }

  export async function asyncValidatorDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:AsyncValidatorDemoFactory', { params });
  }

  export async function globalValidatorDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:GlobalValidatorDemoFactory', { params });
  }

  export async function validatorActionsDemoFactory(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:ValidatorActionsDemoFactory', { params });
  }

  export async function testViewerComponent(): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestViewerComponent', {});
  }

  export async function testFromComponent(): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestFromComponent', {});
  }

  export async function testElements(): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestElements', {});
  }

  export async function testAdd2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestAdd2', { a, b });
  }

  export async function testSub2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestSub2', { a, b });
  }

  export async function testMul2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestMul2', { a, b });
  }

  export async function testDiv2(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestDiv2', { a, b });
  }

  export async function testDF1(df: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestDF1', { df });
  }

  export async function testAdd2Error(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestAdd2Error', { a, b });
  }

  export async function testMultiarg5(a: number, b: number, c: number, d: number, e: number): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:TestMultiarg5', { a, b, c, d, e });
  }

  export async function mockWrapper1(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:MockWrapper1', { params });
  }

  export async function mockWrapper2(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:MockWrapper2', { params });
  }

  export async function mockWrapper3(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:MockWrapper3', { params });
  }

  export async function mockWrapper4(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:MockWrapper4', { params });
  }

  export async function mockWrapper5(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/lib-tests:MockWrapper5', { params });
  }
}
