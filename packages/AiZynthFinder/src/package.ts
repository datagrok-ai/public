/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';

export const _package = new DG.Package();


//top-menu: Chem | RetroSynthesis
//name: RetroSynthesis
//input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
export async function retroSynthesis(molecule: string): Promise<void> {
  const container = await grok.dapi.docker.dockerContainers.filter('aizynthfinder').first();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/aizynthfind', {
    method: 'POST',
    body: JSON.stringify({smiles: molecule}),
    headers: {'Content-Type': 'application/json'},
  });

  console.log((await response.json()));
}