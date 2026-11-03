import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

/** Whether an RDKit-written MOL block is a V3000 one (its counts line says so): RDKit writes V3000 where V2000 cannot
 * hold the molecule (enhanced stereo, more than 999 atoms). */
export function isV3000Block(block: string): boolean {
  let pos = 0;
  for (let k = 0; k < 3 && pos >= 0; k++)
    pos = block.indexOf('\n', pos) + 1;
  const end = block.indexOf('\n', pos);
  return pos > 0 && block.substring(pos, end < 0 ? block.length : end).trimEnd().endsWith('V3000');
}

/** A V3000 MOL block's atom types (MolfileHandler's) and the sum of their `CHG=` charges (elemental analysis:
 * `getAtomsColumn`, the Crux workers' `elementalAnalysis`). */
export function v3000Atoms(block: string): {symbols: string[], charge: number} {
  const handler = MolfileHandler.getInstance(block);
  const charge = handler.getAtomLines()
    .reduce((sum, line) => sum + parseInt(/\sCHG=(-?\d+)/.exec(line)?.[1] ?? '0'), 0);
  return {symbols: handler.atomTypes, charge};
}
