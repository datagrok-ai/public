class ChemPackageDetectors extends DG.Package {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectRDSmiles(col) {
        if (col.name === 'smiles' && col.type === DG.TYPE.STRING) {
            col.semType = DG.SEMTYPE.MOLECULE;
            return col.semType;
        }

        if (col.name === 'rdkit' && col.type === DG.TYPE.STRING) {
            col.semType = 'RDMolecule';
            return col.semType;
        }

        return null;
    }
}
