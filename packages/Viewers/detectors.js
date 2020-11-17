class ViewersPackageDetectors extends DG.Package {

    //tags: semTypeDetector
    //input: column col
    //output: string semType
    detectFlags(col) {
        return (col.type === DG.TYPE.STRING && col.name === 'flag') ? 'flag' : null;
    }
}
