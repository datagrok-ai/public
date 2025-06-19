import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace scripts {

}

export namespace queries {
  export async function dummyPackageQuery(x: number): Promise<DG.DataFrame> {
    return await grok.data.query('APITests:DummyPackageQuery', { x });
  }
}

export namespace funcs {
  export async function getTable(name: string, path: string): Promise<any> {
    return await grok.functions.call('APITests:GetTable', { name, path });
  }

  export async function getColumn(table: DG.DataFrame, columnName: string): Promise<any> {
    return await grok.functions.call('APITests:GetColumn', { table, columnName });
  }

  export async function getDT(rows: number, name: string): Promise<any> {
    return await grok.functions.call('APITests:GetDT', { rows, name });
  }

  export async function getCell(table: DG.DataFrame, rowIndex: number, columnName: string): Promise<any> {
    return await grok.functions.call('APITests:GetCell', { table, rowIndex, columnName });
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('APITests:ExpectTable', { actual, expected });
  }

  export async function dummyPackageFunctionWithDefaultValue(a: string): Promise<any> {
    return await grok.functions.call('APITests:DummyPackageFunctionWithDefaultValue', { a });
  }

  export async function dummyPackageFunction(a: number, b: number): Promise<any> {
    return await grok.functions.call('APITests:DummyPackageFunction', { a, b });
  }

  export async function dummyDataFrameFunction(table: DG.DataFrame): Promise<any> {
    return await grok.functions.call('APITests:DummyDataFrameFunction', { table });
  }

  export async function testIntAsync(a: number): Promise<any> {
    return await grok.functions.call('APITests:TestIntAsync', { a });
  }

  export async function customStringInput(params: any): Promise<any> {
    return await grok.functions.call('APITests:CustomStringInput', { params });
  }

  export async function testOutputAnnotationJoinDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationJoinDf', { data, col });
  }

  export async function testOutputAnnotationJoinCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationJoinCol', { data, col });
  }

  export async function testOutputAnnotationJoinColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationJoinColList', { data, col });
  }

  export async function testOutputAnnotationReplaceDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationReplaceDf', { data, col });
  }

  export async function testOutputAnnotationReplaceCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationReplaceCol', { data, col });
  }

  export async function testOutputAnnotationReplaceColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputAnnotationReplaceColList', { data, col });
  }

  export async function testOutputWithoutAction(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestOutputWithoutAction', { data, col });
  }

  export async function expectDate(actual: any, expected: any): Promise<any> {
    return await grok.functions.call('APITests:ExpectDate', { actual, expected });
  }

  export async function testVectorFunc(col: DG.Column, prefix: string): Promise<any> {
    return await grok.functions.call('APITests:TestVectorFunc', { col, prefix });
  }

  export async function testVectorFuncNonVectorizableParam(col: DG.Column, prefix: string, postfix: DG.Column): Promise<any> {
    return await grok.functions.call('APITests:TestVectorFuncNonVectorizableParam', { col, prefix, postfix });
  }
}
