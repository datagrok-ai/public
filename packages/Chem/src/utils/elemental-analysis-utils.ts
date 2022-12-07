import * as DG from 'datagrok-api/dg';
import { V2000_ATOM_NAME_LEN, V2000_ATOM_NAME_POS } from '../constants';
import { convertMolNotation } from '../package';

/** Gets map of chem elements to list with counts of atoms in rows */
export function getAtomsColumn(molCol: DG.Column): [Map<string, Int32Array>, number[]] {
    let elements: Map<string, Int32Array> = new Map();
    const invalid: number[] = new Array<number>();
    let smiles = molCol.getTag(DG.TAGS.UNITS) === DG.UNITS.Molecule.SMILES;
    for (let rowI = 0; rowI < molCol.length; rowI++) {
      let el: string = molCol.get(rowI);
      if (smiles) {
        try {
          el = convertMolNotation(el, DG.UNITS.Molecule.SMILES, DG.UNITS.Molecule.MOLBLOCK);
        } 
        catch {
          invalid.push(rowI);
        }
      }
      let curPos = 0;
      curPos = el.indexOf('\n', curPos) + 1;
      curPos = el.indexOf('\n', curPos) + 1;
      curPos = el.indexOf('\n', curPos) + 1;

      const atomCounts = parseInt(el.substring(curPos, curPos + 3));
  
      for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
        curPos = el.indexOf('\n', curPos) + 1;
        const elName: string = el
          .substring(curPos + V2000_ATOM_NAME_POS, curPos + V2000_ATOM_NAME_POS + V2000_ATOM_NAME_LEN)
          .trim();
  
        if (!elements.has(elName))
          elements.set(elName, new Int32Array(molCol.length));
          
        ++elements.get(elName)![rowI];
      } 
    }
    return [elements, invalid];
  }

/** Check that packages are installed */
export function checkPackage(packageName: string, functionName: string) : boolean {
  const funcList: DG.Func[] = DG.Func.find({package: packageName, name: functionName});
  console.debug(`${packageName}: ${functionName} funcList.length = ${funcList.length}`);
  return funcList.length === 1 ? true : false;
}

