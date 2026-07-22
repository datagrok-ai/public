/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.func()
  static info() : void {
    grok.shell.info(_package.webRoot);
  }
}

//name: getColumn
//input: dataframe table
//input: string columnName
//output: column col
export function getColumn(table: DG.DataFrame, columnName: string): DG.Column {
  const col = table.getCol(columnName);
  return col;
}

// ---- JS queue functions: run server-side in the celery Node worker (meta.queue) ----

//name: jsCvmInt
//meta.queue: true
//input: int x
//output: int result
export function jsCvmInt(x: number): number {
  return x;
}

//name: jsCvmDouble
//meta.queue: true
//input: double x
//output: double result
export function jsCvmDouble(x: number): number {
  return x;
}

//name: jsCvmBool
//meta.queue: true
//input: bool x
//output: bool result
export function jsCvmBool(x: boolean): boolean {
  return x;
}

//name: jsCvmString
//meta.queue: true
//input: string x
//output: string result
export function jsCvmString(x: string): string {
  return x;
}

//name: jsCvmBigInt
//meta.queue: true
//input: bigint x
//output: bigint result
export function jsCvmBigInt(x: bigint): bigint {
  return x;
}

//name: jsCvmDataframe
//meta.queue: true
//input: dataframe df
//output: dataframe result
export function jsCvmDataframe(df: DG.DataFrame): DG.DataFrame {
  return df;
}

//name: jsCvmEmptyDataframe
//meta.queue: true
//output: dataframe result
export function jsCvmEmptyDataframe(): DG.DataFrame {
  return DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'col1', [])]);
}

//name: jsCvmError
//meta.queue: true
//output: string result
export function jsCvmError(): string {
  throw new Error('planned jsCvmError failure');
}

//name: jsCvmCancel
//meta.queue: true
//input: int seconds
//output: string result
export async function jsCvmCancel(seconds: number): Promise<string> {
  await new Promise((resolve) => setTimeout(resolve, seconds * 1000));
  return 'done';
}

//name: jsCvmProgress
//meta.queue: true
//output: string result
export function jsCvmProgress(): string {
  (globalThis as any).DG_TASK_PROGRESS?.(50, 'half way');
  return 'ok';
}

//name: jsCvmCurrentUser
//meta.queue: true
//output: string result
export async function jsCvmCurrentUser(): Promise<string> {
  return (await grok.dapi.users.current()).login;
}

//name: jsCvmServer
//meta.server: true
//input: string x
//output: string result
export function jsCvmServer(x: string): string {
  return x;
}

//name: jsCvmCustomContainer
//meta.queue: node-worker
//output: string result
export function jsCvmCustomContainer(): string {
  return 'custom';
}
