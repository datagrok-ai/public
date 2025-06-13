import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace queries {
  export async function dummyPackageQuery(x: number): Promise<DG.DataFrame> {
    return await grok.data.query('@datagrok/api-tests:DummyPackageQuery', { x });
  }
}

export namespace funcs {
  export async function getTable(name: string, path: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:GetTable', { name, path });
  }

  export async function getColumn(table: DG.DataFrame, columnName: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:GetColumn', { table, columnName });
  }

  export async function getDT(rows: number, name: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:GetDT', { rows, name });
  }

  export async function getCell(table: DG.DataFrame, rowIndex: number, columnName: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:GetCell', { table, rowIndex, columnName });
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:ExpectTable', { actual, expected });
  }

  export async function dummyPackageFunctionWithDefaultValue(a: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:DummyPackageFunctionWithDefaultValue', { a });
  }

  export async function dummyPackageFunction(a: number, b: number): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:DummyPackageFunction', { a, b });
  }

  export async function dummyDataFrameFunction(table: DG.DataFrame): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:DummyDataFrameFunction', { table });
  }

  export async function testIntAsync(a: number): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestIntAsync', { a });
  }

  export async function customStringInput(params: any): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:CustomStringInput', { params });
  }

  export async function testOutputAnnotationJoinDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationJoinDf', { data, col });
  }

  export async function testOutputAnnotationJoinCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationJoinCol', { data, col });
  }

  export async function testOutputAnnotationJoinColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationJoinColList', { data, col });
  }

  export async function testOutputAnnotationReplaceDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationReplaceDf', { data, col });
  }

  export async function testOutputAnnotationReplaceCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationReplaceCol', { data, col });
  }

  export async function testOutputAnnotationReplaceColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputAnnotationReplaceColList', { data, col });
  }

  export async function testOutputWithoutAction(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestOutputWithoutAction', { data, col });
  }

  export async function expectDate(actual: any, expected: any): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:ExpectDate', { actual, expected });
  }

  export async function testVectorFunc(col: DG.Column, prefix: string): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestVectorFunc', { col, prefix });
  }

  export async function testVectorFuncNonVectorizableParam(col: DG.Column, prefix: string, postfix: DG.Column): Promise<any> {
    return await grok.functions.call('@datagrok/api-tests:TestVectorFuncNonVectorizableParam', { col, prefix, postfix });
  }
}
