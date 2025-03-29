/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/aizynthfinder.css';
import {AiZynthFinderViewer} from './aizynthfinder-viewer';
import {createPathsTreeTabs, isFragment} from './utils';
import {ReactionData, Tree} from './aizynth-api';

export const _package = new DG.Package();


//name: CalculateRetroSynthesisPaths
//meta.cache: all
//meta.cache.invalidateOn: 0 0 * * *
//input: string molecule = "O=C1Nc2ccccc2C(C2CCCCC2)=NC1" { semType: Molecule }
//output: string paths
export async function calculateRetroSynthesisPaths(molecule: string): Promise<string> {
  const container = await grok.dapi.docker.dockerContainers.filter('retrosynthesis').first();
  const startTime = performance.now();
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(container.id, '/aizynthfind', {
    method: 'POST',
    body: JSON.stringify({smiles: molecule}),
    headers: {'Content-Type': 'application/json'},
  });
  console.log(`Request to aizynthfinder finished in ${performance.now() - startTime} ms`);
  const resJson = await response.json();
  if (resJson[1] !== 200)
    throw new Error('Error occured during paths generation');

  return resJson[0].result;
}


//name: Chemistry | Retrosynthesis
//tags: panel, chem, widgets
//condition: true
//input: string smiles { semType: Molecule }
//output: widget result
export async function retroSynthesisPath(molecule: string): Promise<DG.Widget> {
  if (!molecule || DG.chem.Sketcher.isEmptyMolfile(molecule))
    return new DG.Widget(ui.divText('Molecule is empty'));
  if (DG.chem.isSmarts(molecule) || isFragment(molecule))
    return new DG.Widget(ui.divText('Not applicable for smarts or moleculer fragments'));

  //check molecule is valid and convert to smiles
  try {
    molecule = DG.chem.convert(molecule, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
  } catch {
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
  }

  // Call retrosynthesis function
  let result: string;
  try {
    result = await grok.functions.call('Retrosynthesis:calculateRetroSynthesisPaths',
      {molecule: molecule});
  } catch (e: any) {
    return new DG.Widget(ui.divText(e));
  }

  // Parse and process reaction data
  try {
    const reactionData: ReactionData = JSON.parse(result);
    const paths: Tree[] = reactionData?.data?.[0]?.trees;
    return paths?.length ?
      new DG.Widget(createPathsTreeTabs(paths, false).root) :
      new DG.Widget(ui.divText('No paths found for the molecule'));
  } catch {
    return new DG.Widget(ui.divText('Error processing retrosynthesis data'));
  }
}


//name: Retrosynthesis Viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/chem-similarity-search-viewer.svg
export function retrosynthesisViewer(): AiZynthFinderViewer {
  return new AiZynthFinderViewer();
}

//top-menu: Chem | Retrosynthesis
//name: RetrosynthesisPath
export function retrosynthesisTopMenu(): void {
  (grok.shell.v as DG.TableView).addViewer('Retrosynthesis Viewer');
}

