class SequencetranslatorPackageDetectors extends DG.Package {
  // https://stackoverflow.com/questions/19269545/how-to-get-a-number-of-random-elements-from-an-array
  getRandom(arr, n) {
    let result = new Array(n),
      len = arr.length,
      taken = new Array(len);
    while (n--) {
      let x = Math.floor(Math.random() * len);
      result[n] = arr[x in taken ? taken[x] : x];
      taken[x] = --len in taken ? taken[len] : len;
    }
    return result;
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      const numberOfRandomValues = Math.max(1, Math.round(0.1 * col.categories.length));
      const decimatedColumn = this.getRandom(col.categories.slice(0, col.categories.length - 1), numberOfRandomValues);
      if (this.isClassicFormat(decimatedColumn)) return 'dna_sequence/classic code';
      if (this.isBioSpringCodeForASOGapmers(decimatedColumn)) return 'BioSpring Code For ASO Gapmers';
      if (this.isGCRSCodeForASOGampers(decimatedColumn)) return 'Janssen GCRS code For ASO Gapmers';
    }
  }

  isClassicFormat(decimatedColumn) {
    for (let valueInRow of decimatedColumn)
      if (!/^[ATGC]+$/.test(valueInRow)) return false;
    return true;
  }

  isBioSpringCodeForASOGapmers(decimatedColumn) {
    for (let valueInRow of decimatedColumn)
      if (!/^[*56789ACTG]+$/.test(valueInRow)) return false;
    return true;
  }

  isGCRSCodeForASOGampers(decimatedColumn) {
    for (let valueInRow of decimatedColumn)
      if (!(valueInRow.slice(0, 3) === 'moe')) return false;
    return true;
  }
}