import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package';

export async function nglViewerGen(): Promise<void> {
  const mol2Str: string = await _package.files.readAsText('samples/dc.mol2');

  const atomCount = 76;
  const x: number[] = new Array<number>(atomCount).fill(0);
  const y: number[] = new Array<number>(atomCount).fill(0);
  const z: number[] = new Array<number>(atomCount).fill(0);

  const resMol2Str = replaceCoordsToMol2(mol2Str, x, y, z);

  console.debug(resMol2Str);
}

const TRIPOS_ATOM_LINE = '@<TRIPOS>ATOM';
const TRIPOS_BOND_LINE = '@<TRIPOS>BOND';

function replaceCoordsToMol2(mol2: string, x: number[], y: number[], z: number[]): string {
  const atomLinePos = mol2.indexOf(TRIPOS_ATOM_LINE);
  if (atomLinePos == -1) throw new Error('Invalid data');
  const bondLinePos = mol2.indexOf(TRIPOS_BOND_LINE, atomLinePos);
  if (bondLinePos == -1) throw new Error('Invalid data');

  const res: string = '';
  return res;
}