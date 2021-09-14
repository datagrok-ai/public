import * as DG from 'datagrok-api/dg';

class PeptidesPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col: DG.Column) {

    if (col.type === DG.TYPE.STRING) {
      for (let i = 0; i != col.categories.length; i++) {
        let papRgx = new RegExp("^((\\w+)?-){5,49}\\w+$").test(col.categories[i]);
        if (!papRgx)
          return null;
      }
      col.semType = 'alignedSequence';
      return (col.semType);
    }

    return null;
  }
}
