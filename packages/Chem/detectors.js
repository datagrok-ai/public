class ChemsPackageDetectors extends DG.Package {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectSmiles(col) {
        if (col.name === 'smiles' && col.type === DG.TYPE.STRING) {
            col.semType = DG.SEMTYPE.MOLECULE;
            return col.semType;
        }

        return null;
    }
}
