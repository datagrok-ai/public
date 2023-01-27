import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface IPdbHelper {
  pdbToDf(pdb: string, name: string): Promise<DG.DataFrame>;
}

export async function getPdbHelper(): Promise<IPdbHelper> {
  const packageName = 'BiostructureViewer';
  const funcList = DG.Func.find({package: packageName, name: 'getPdbHelper'});
  if (funcList.length === 0)
    throw new Error(`Package '${packageName}' must be installed for PdbHelper.`);

  const res: IPdbHelper = (await funcList[0].prepare().call()).getOutputParamValue() as IPdbHelper;
  return res;
}