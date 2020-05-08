class ChemspacePackageDetectors extends GrokPackage {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectCsId(col) {
        if (col.name === 'CS-id') {
            col.semType = 'chemspace-id';
            return col.semType;
        }

        return null;
    }
}
