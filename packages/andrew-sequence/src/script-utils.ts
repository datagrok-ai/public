import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export type ScriptType = 'JS' | 'Python';

export async function selectOredersByCountry(packageName: string, country: string): Promise<DG.DataFrame> {
  return await grok.data.query(`${packageName}:ordersByCountry`, {country});
}

export async function countSubsequences(
  packageName: string, params: { [key: string]: unknown }, scriptType: ScriptType = 'Python',
): Promise<number> {
  const scriptName = `${packageName}:CountSubsequence${scriptType}`;
  return await grok.functions.call(scriptName, params);
}

export async function getSubsequencesCountColumn(
  packageName: string, params: { [key: string]: unknown }, scriptType: ScriptType = 'Python',
): Promise<DG.Column> {
  const scriptName = `${packageName}:CountSubsequence${scriptType}Dataframe`;
  const df = await grok.functions.call(scriptName, params);
  const countCol = df.columns.byIndex(0);
  return countCol;
}
