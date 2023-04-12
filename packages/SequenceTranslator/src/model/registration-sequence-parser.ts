const DUPLEX_PREFIX = {
  SS: 'SS ',
  AS: '\nAS ',
};

const DIMER_PREFIX = {
  SS: 'SS ',
  AS1: '\nAS1 ',
  AS2: '\nAS2 ',
};

export class RegistrationSequenceParser {
  getDuplexStrands(strands: string): {ss: string, as: string} {
    const result = this.prepareInput(strands).slice(DUPLEX_PREFIX.SS.length).split(DUPLEX_PREFIX.AS);
    return {ss: result[0], as: result[1]};
  }

  getDimerStrands(strands: string): {ss: string, as1: string, as2: string} {
    let result = this.prepareInput(strands).slice(DIMER_PREFIX.SS.length).split(DIMER_PREFIX.AS1);
    const as = result[0];
    result = result[1].split(DIMER_PREFIX.AS2);
    const as1 = result[0];
    const as2 = result[1];
    return {ss: as, as1: as1, as2: as2};
  }

  // todo: port validation to dedicated class
  private prepareInput(strands: string) {
    return strands.replace(/\r/g, '');
  }
}
