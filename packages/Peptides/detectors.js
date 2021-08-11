class PeptidesPackageDetectors extends DG.Package {

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {

    if (col.type === DG.TYPE.STRING) {
      for (let i = 0; i != col.categories.length; i++) {
        var pepRgx = new RegExp(
          "(G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T|(\[(G|L|Y|S|E|Q|D|N|F|A|K|R|H|C|V|P|W|I|M|T).+\])){3,30}")
          .test(col.categories[i]);
        if (!pepRgx)
          return null;
      }
      col.semType = 'peptide';
      return (col.semType);
    }

    return null;

  }
}