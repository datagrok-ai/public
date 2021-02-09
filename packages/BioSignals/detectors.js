class BiosignalsPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectBioS(col) {
    if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) && (col.name.match(/ecg/i) ||
      col.name.match(/eeg/i) || col.name.match(/eda/i))) {
      col.semType = 'BioSignal-' + col.name.toLowerCase();
      return col.semType;
    }

    return null;
  }

  //input: dataframe table
  //output: bool result
  analysisCondition(table) {
    let columns = table.columns.toList();

    // for (let column in columns) {
    //   grok.shell.info(column.semType.includes('BioSignal'));
    // }

    //return columns.some((c) => c.semType.includes('BioSignal'));
    return true;
  }
}
