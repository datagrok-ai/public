import {_rdKitModule} from "./chem-common-rdkit";

function syncQueryAromatics_1(molBlock: string,  bonds2Change: Array<number> | null = null) : string | Array<number> {
  let curPos = 0;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  curPos = molBlock.indexOf('\n', curPos) + 1;
  const atomCounts = parseInt(molBlock.substring(curPos, curPos + 3));
  const bondCounts = parseInt(molBlock.substring(curPos + 3, curPos + 6));

  for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
  }

  const read = bonds2Change === null;
  bonds2Change ??= [];

  let bondOrder = -1;
  for (let bondRowI = 0; bondRowI < bondCounts; bondRowI++) {
    curPos = molBlock.indexOf('\n', curPos) + 1;
    if (read) {
      bondOrder = parseInt(molBlock.substring(curPos + 8, curPos + 9));
      if (bondOrder === 4)
        bonds2Change.push(bondRowI);
    }
    else {
      if (bonds2Change.includes(bondRowI))
        molBlock = molBlock.slice(0, curPos + 8) + '4' + molBlock.slice(curPos + 9);
    }
  }
  return read ? bonds2Change : molBlock;
}

function syncQueryAromatics_2(molBlockAroma: string, molBlock : string) : string {
  const bonds2Change = syncQueryAromatics_1(molBlock);
  const molModified = syncQueryAromatics_1(molBlockAroma, bonds2Change as Array<number>);
  return molModified as string;
}

export function aromatizeMolBlock(molString: string) : string {
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
  const newQueryMolString = syncQueryAromatics_2(molBlockAroma, molString);
  return newQueryMolString;
}