class BiosignalsPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectBioS(col) {
    if (col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) {
      const columnName = col.name.toLowerCase();
      if (columnName.match(/ecg/i)) {
        col.semType = 'BioSignal-ECG';
      } else if (columnName.match(/eda/i)) {
        col.semType = 'BioSignal-EDA';
      } else if (columnName.match(/accel/i)) {
        col.semType = 'BioSignal-Accelerometer';
      } else if (columnName.match(/emg/i)) {
        col.semType = 'BioSignal-EMG';
      } else if (columnName.match(/eeg/i)) {
        col.semType = 'BioSignal-EEG';
      } else if (columnName.match(/abp/i)) {
        col.semType = 'BioSignal-ABP';
      } else if (columnName.match(/bvp/i) || columnName.match(/ppg/i)) {
        col.semType = 'BioSignal-BVP(PPG)';
      } else if (columnName.match(/resp/i)) {
        col.semType = 'BioSignal-Respiration';
      }
      return col.semType;
    }
    return null;
  }

  //input: dataframe table
  //output: bool result
  analysisCondition(table) {
    let columns = table.columns.toList();

    function check(semType) {
      return columns.some((c) => c.semType === semType)
    }

    return (check('BioSignal-ECG') || check('BioSignal-EDA') || check('BioSignal-Accelerometer')
      || check('BioSignal-EMG') || check('BioSignal-EEG') || check('BioSignal-ABP') ||
      check('BioSignal-BVP(PPG)') || check('BioSignal-Respiration'));
  }
}
