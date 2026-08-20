import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';
import {V2K_CONST} from '@datagrok-libraries/chem-meta/src/formats/molfile-const';

const ATOM_LINE_TAIL = ' 0  0  0  0  0  0  0  0  0  0  0  0';

const pad3 = (v: number): string => String(v).padStart(3);
const padCoord = (v: number): string => v.toFixed(4).padStart(10);

/** Rebuilds a V3000 molfile as V2000 (coordinates, elements, bonds) for parsers lacking V3000 support (Mol* 2.x). */
export function molfileV3KToV2K(molfileV3K: string): string {
  const handler = MolfileHandler.getInstance(molfileV3K);
  const [atomCount, bondCount] = [handler.atomCount, handler.bondCount];
  if (atomCount > V2K_CONST.MAX_ATOM_COUNT || bondCount > V2K_CONST.MAX_ATOM_COUNT) {
    throw new Error(
      `V3000 molfile with over ${V2K_CONST.MAX_ATOM_COUNT} atoms or bonds is not convertible to V2000.`);
  }

  const srcLines = molfileV3K.replace(/\r/g, '').split('\n');
  const lines: string[] = [srcLines[0] ?? '', srcLines[1] ?? '', srcLines[2] ?? ''];
  lines.push(`${pad3(atomCount)}${pad3(bondCount)}  0  0  0  0  0  0  0  0999 V2000`);
  const [x, y, z, atomTypes] = [handler.x, handler.y, handler.z, handler.atomTypes];
  for (let i = 0; i < atomCount; i++)
    lines.push(`${padCoord(x[i])}${padCoord(y[i])}${padCoord(z[i])} ${atomTypes[i].padEnd(3)}${ATOM_LINE_TAIL}`);
  const [pairs, bondTypes] = [handler.pairsOfBondedAtoms, handler.bondTypes];
  for (let i = 0; i < bondCount; i++)
    lines.push(`${pad3(pairs[i][0])}${pad3(pairs[i][1])}${pad3(bondTypes[i])}  0`);
  lines.push(V2K_CONST.END);
  return lines.join('\n') + '\n';
}
