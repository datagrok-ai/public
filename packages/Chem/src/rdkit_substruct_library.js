class RdKitSubstructLibrary {

    constructor(module) {

        this.rdKitModule = module;
        this.library = null;

    }

    init(dict) {

        this.deinit();
        this.library = new this.rdKitModule.SubstructLibrary();
        for (let item of dict) {
            let mol = null;
            try {
                mol = this.rdKitModule.get_mol(item);
                this.library.add_mol(mol);
                mol.delete();
            } catch (e) {
                console.error(
                    "Possibly a malformed molString: `" + item + "`");
                // Won't rethrow
            }
        }

    }

    search(query) {

        const queryMol = this.rdKitModule.get_mol(query);
        const matches = this.library.get_matches(queryMol, false, 1, -1);
        queryMol.delete();
        return matches;

    }

    deinit() {

        this.library?.delete();

    }

}