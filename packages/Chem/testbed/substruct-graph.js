function searchByGraph(inputSmilesSet, querySmiles) {
  
  let molQuery = Module.get_mol(querySmiles);
  for (let smilesSet of inputSmilesSet) {
    let molSet = Module.get_mol(smilesSet);
    let match = molSet.get_substruct_match(molQuery);
    molSet.delete();
  }
  molQuery.delete();
  
}