import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {

}

export namespace Queries {
  export async function dummyPackageQuery(x: number): Promise<DG.DataFrame> {
    return await grok.data.query('ApiTests:DummyPackageQuery', { x });
  }
}

export namespace Funcs {
  export async function getTable(name: string, path: string): Promise<any> {
    return await grok.functions.call('ApiTests:GetTable', { name, path });
  }

  export async function getColumn(table: DG.DataFrame, columnName: string): Promise<any> {
    return await grok.functions.call('ApiTests:GetColumn', { table, columnName });
  }

  export async function getDT(rows: number, name: string): Promise<any> {
    return await grok.functions.call('ApiTests:GetDT', { rows, name });
  }

  export async function getCell(table: DG.DataFrame, rowIndex: number, columnName: string): Promise<any> {
    return await grok.functions.call('ApiTests:GetCell', { table, rowIndex, columnName });
  }

  export async function expectTable(actual: DG.DataFrame, expected: DG.DataFrame): Promise<any> {
    return await grok.functions.call('ApiTests:ExpectTable', { actual, expected });
  }

  export async function dummyPackageFunctionWithDefaultValue(a: string): Promise<any> {
    return await grok.functions.call('ApiTests:DummyPackageFunctionWithDefaultValue', { a });
  }

  export async function dummyPackageFunction(a: number, b: number): Promise<any> {
    return await grok.functions.call('ApiTests:DummyPackageFunction', { a, b });
  }

  export async function dummyDataFrameFunction(table: DG.DataFrame): Promise<any> {
    return await grok.functions.call('ApiTests:DummyDataFrameFunction', { table });
  }

  export async function testIntAsync(a: number): Promise<any> {
    return await grok.functions.call('ApiTests:TestIntAsync', { a });
  }

  export async function customStringInput(params: any): Promise<any> {
    return await grok.functions.call('ApiTests:CustomStringInput', { params });
  }

  export async function testOutputAnnotationJoinDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationJoinDf', { data, col });
  }

  export async function testOutputAnnotationJoinCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationJoinCol', { data, col });
  }

  export async function testOutputAnnotationJoinColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationJoinColList', { data, col });
  }

  export async function testOutputAnnotationReplaceDf(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationReplaceDf', { data, col });
  }

  export async function testOutputAnnotationReplaceCol(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationReplaceCol', { data, col });
  }

  export async function testOutputAnnotationReplaceColList(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputAnnotationReplaceColList', { data, col });
  }

  export async function testOutputWithoutAction(data: DG.DataFrame, col: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestOutputWithoutAction', { data, col });
  }

  export async function expectDate(actual: any, expected: any): Promise<any> {
    return await grok.functions.call('ApiTests:ExpectDate', { actual, expected });
  }

  export async function testVectorFunc(col: DG.Column, prefix: string): Promise<any> {
    return await grok.functions.call('ApiTests:TestVectorFunc', { col, prefix });
  }

  export async function testVectorFuncNonVectorizableParam(col: DG.Column, prefix: string, postfix: DG.Column): Promise<any> {
    return await grok.functions.call('ApiTests:TestVectorFuncNonVectorizableParam', { col, prefix, postfix });
  }
}
