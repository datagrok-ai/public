class RdKitSubstructLibrary {

  constructor(module) {

    this.rdKitModule = module;
    this.library = null;

  }

  init(dict) {

    this.deinit();
    if (dict.length === 0) { this.library = null; return; }
    this.library = new this.rdKitModule.SubstructLibrary();
    for (let item of dict) {
      try {
        let mol = this.rdKitModule.get_mol(item);
        this.library.add_mol(mol);
        mol.delete();
      } catch (e) {
        console.error(
          "Possibly a malformed molString: `" + item + "`");
        // preserving indices with a placeholder
        let mol = this.rdKitModule.get_mol('');
        this.library.add_mol(mol);
        mol.delete();
        // Won't rethrow
      }
    }

  }

  search(query) {

    if (this.library == null) { return "[]"; }
    const queryMol = this.rdKitModule.get_mol(query);
    queryMol.merge_hs_as_queries();
    const matches = this.library.get_matches(queryMol, false, 1, -1);
    queryMol.delete();
    return matches;

  }

  deinit() {

    this.library?.delete();
    this.library = null;

  }

}