import {RDModule} from "@datagrok-libraries/chem-meta/src/rdkit-api";

function fillAromaticBondsToChange(molBlock: string) : number[] {
  let curPos = 0;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  const atomCounts = parseInt(molBlock.substring(curPos, curPos + 3));
  const bondCounts = parseInt(molBlock.substring(curPos + 3, curPos + 6));

  for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
  }

  const bonds2Change : number[] = [];
  let bondOrder = -1;
  for (let bondRowI = 0; bondRowI < bondCounts; bondRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
    bondOrder = parseInt(molBlock.substring(curPos + 8, curPos + 9));
    if (bondOrder === 4)
      bonds2Change.push(bondRowI);
  }
  return bonds2Change;
}

function changeToAromaticsBonds(molBlock: string,  bonds2Change: Array<number>) : string {
  let curPos = 0;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  const atomCounts = parseInt(molBlock.substring(curPos, curPos + 3));
  const bondCounts = parseInt(molBlock.substring(curPos + 3, curPos + 6));

  for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
  }

  for (let bondRowI = 0; bondRowI < bondCounts; bondRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
    if (bonds2Change.includes(bondRowI))
      molBlock = molBlock.slice(0, curPos + 8) + '4' + molBlock.slice(curPos + 9);
  }
  return molBlock;
}

export function syncQueryAromatics(molBlockAroma: string, molBlock : string) : string {
  const bonds2Change = fillAromaticBondsToChange(molBlock);
  return changeToAromaticsBonds(molBlockAroma, bonds2Change);
}

export function aromatizeMolBlock(molString: string, _rdKitModule: RDModule) : string {
  let molTmp = null;
  try { molTmp = _rdKitModule.get_mol(molString); }
  catch(e) {
    if (molTmp !== null && molTmp.is_valid())
      molTmp.delete();

    try { molTmp = _rdKitModule.get_qmol(molString); }
    catch(e) {
      return molString;
    }
  }
  let molBlockAroma = null;
  try { molBlockAroma = molTmp!.get_aromatic_form(); }
  catch(e) {
    molBlockAroma = molString;
  }

  molTmp!.delete();
  const newQueryMolString = syncQueryAromatics(molBlockAroma, molString);
  return newQueryMolString;
}