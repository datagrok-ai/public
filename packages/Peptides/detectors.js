class PeptidesPackageDetectors extends DG.Package {
    
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectAligned(col) {

    if (col.type === DG.TYPE.STRING) {
      for (let i = 0; i != col.categories.length; i++) {
        //var papRgx = new RegExp("^((\w+)?-){5,49}\w+$");

        var papRgx = new RegExp("NH2").test(col.categories[i]);

        if (!papRgx)
          return null;
      }
      col.semType = 'alignedSequence';
      return (col.semType);
    }
  }
}
