class ChartsPackageDetectors extends DG.Package {
  //meta.role: semTypeDetector
  //input: column col {type: string}
  //output: string semType
  detectFlags(col) {
    return (col.type === DG.TYPE.STRING && col.name === 'flag') ? 'flag' : null;
  }

  //meta.role: semTypeDetector
  //input: column col {type: numerical}
  //output: string semType
  detectMagnitude(col) {
    if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) &&
      (0 < col.min && col.max < 10) && col.name.toLowerCase() === 'magnitude')
      return 'Magnitude';
    return null;
  }
}
