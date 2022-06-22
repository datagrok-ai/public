var library = null;

function buildLibrary(smilesArray)
{
  if (library != null) {
    library.delete();
  }
  library = new Module.SubstructLibrary();
  for (smiles of smilesArray) {
      try {
        library.add_smiles(smiles);
      } catch (e) {
      }
  }

}

function searchInLibrary(smilesString)
{
  let query = Module.get_mol(smilesString);
  const matches = library.get_matches(query);
  query.delete();
}

function closeLibrary() {
  library.delete();
  library = null;
}