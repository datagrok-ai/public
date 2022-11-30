import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export async function findNewick(data: DG.DataFrame): Promise<string> {
  return await grok.functions.call('PhyloTreeViewer:getNewick', {data});
}