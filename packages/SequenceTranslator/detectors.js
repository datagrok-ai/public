class SequencetranslatorPackageDetectors extends DG.Package {
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectNucleotides(col) {
    if (col.type === DG.TYPE.STRING) {
      if (this.isClassicFormat(col)) {
        return 'dna_sequence/classic';
      } else if (this.isRobotFormat(col)) {
        col.semType = 'dna_sequence/robot';
        return col.semType;
      } else if (this.isCROFormat(col)) {
        col.semType = 'dna_sequence/cro';
        return col.semType;
      }
    }
    return null;
  }

  isClassicFormat(col) {
    for (let valueInRow of col.categories) {
      if (!this.checkWhetherClassicFormat(valueInRow)){ return false; }
    }
    return true;
  }

  //input: string nucleotide
  //output: bool result
  checkWhetherClassicFormat(inputString) {
    const dnaBases = new Set(['A','C','T','G']);
    return inputString.split('').every(x => dnaBases.has(x));
  }

  isRobotFormat(col) {
    for (let valueInRow of col.categories) {
      if (!this.checkWhetherRobotFormat(valueInRow)){ return false; }
    }
    return true;
  }

  //input: string nucleotide
  //output: bool result
  checkWhetherRobotFormat(inputString) {
    const dnaBases = new Set(['*','5','6','7','8','9','A','C','T','G']);
    return inputString.split('').every(x => dnaBases.has(x));
  }

  isCROFormat(col) {
    for (let valueInRow of col.categories) {
      if (!this.checkWhetherCROFormat(valueInRow)){ return false; }
    }
    return true;
  }

  //input: string nucleotide
  //output: bool result
  checkWhetherCROFormat(inputString) {
    return inputString.slice(0, 3) === 'moe';
  }
}
