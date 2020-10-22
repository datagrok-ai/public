class BiosensorsPackageDetectors extends DG.Package {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectBioS(col) {
        if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) && (col.name.match(/ecg/i) ||
            col.name.match(/eeg/i) || col.name.match(/eda/i))) {
            col.semType = 'Biosensor-' + col.name.toLowerCase();
            return col.semType;
        }

        return null;
    }

    //input: dataframe table
    //output: bool result
    analysisCondition(table) {
        let columns = table.columns.toList();

        function checkSemType(semType) {
            return columns.some((c) => c.semType.includes(semType))
        }

        if (checkSemType('Biosensor')) {
            grok.shell.info('check successful')
            return true;
        }
    }

}
