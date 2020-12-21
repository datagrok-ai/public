class ViewersPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectFlags(col) {
    return (col.type === DG.TYPE.STRING && col.name === 'flag') ? 'flag' : null;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMagnitude(col) {
    if ((col.type === DG.TYPE.FLOAT || col.type === DG.TYPE.INT) &&
      (0 < col.min && col.max < 10) && col.name.toLowerCase() === 'magnitude') {
      col.semType = 'Magnitude';
      return col.semType;
    }
    return null;
  }
}
