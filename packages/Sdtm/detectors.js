class SDTMPackageDetectors extends DG.Package {
    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detector(col) {
        if (col.dataFrame.name === 'stdmlb' && col.name.toLowerCase() === 'studyid')
            return col.semType = 'Study';

        if (col.dataFrame.name === 'stdmlb' && col.name.toLowerCase() === 'usubjid')
            return col.semType = 'Subject';

        return null;
    }
}
