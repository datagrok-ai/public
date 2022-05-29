function searchByGraph(inputSmilesSet, querySmiles) {
  let cntr = 0;
  let molQuery = RDKitModule.get_qmol(querySmiles);
  for (let idx = 0; idx < inputSmilesSet.length; ++idx) {
    const molString = inputSmilesSet[idx];
    try {
      let molSet = RDKitModule.get_mol(molString);
      let match = molSet.get_substruct_match(molQuery);
      if (match !== '{}') {
        cntr++;
      }
    molSet.delete();
    } catch (e) {
      console.log(`Cannot process a pattern ${idx}: \r\n`, molString);
    }
  }
  molQuery.delete();
  return cntr;
}