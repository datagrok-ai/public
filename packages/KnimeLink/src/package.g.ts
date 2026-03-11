import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: KNIME
//output: view result
//meta.role: app
//meta.browsePath: Compute
export async function knimeLinkApp() : Promise<any> {
  return await PackageFunctions.knimeLinkApp();
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: KNIME
export async function knimeLinkAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.knimeLinkAppTreeBrowser(treeNode);
}

//name: Execute KNIME Workflow
//input: string workflowId 
//input: string inputJson { nullable: true }
//input: dataframe inputTable { nullable: true }
//input: string tableParamName { nullable: true }
//output: dataframe result
export async function executeKnimeWorkflow(workflowId: string, inputJson?: string, inputTable?: DG.DataFrame, tableParamName?: string) : Promise<any> {
  return await PackageFunctions.executeKnimeWorkflow(workflowId, inputJson, inputTable, tableParamName);
}
