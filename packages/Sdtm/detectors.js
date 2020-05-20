class SDTMPackageDetectors extends grok.Package {
    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detector(col) {
        if (col.name.toLowerCase() === 'studyid')
            return col.semType = 'Study';

        if (col.name.toLowerCase() === 'usubjid')
            return col.semType = 'Subject';

        return null;
    }
}
