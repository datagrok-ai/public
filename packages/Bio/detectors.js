/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 *
 * TODO: Use detectors from WebLogo pickUp.. methods
 */


class BioPackageDetectors extends DG.Package {

  static PeptideFastaAlphabet = new Set([
    'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
    'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
  ]);

  static DnaFastaAlphabet = new Set(['A', 'C', 'G', 'T']);

  static RnaFastaAlphabet = new Set(['A', 'C', 'G', 'U']);

  static SmilesRawAlphabet = new Set([
    'A', 'B', 'C', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'R', 'S', 'Z',
    'a', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'r', 's', 't', 'u',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '+', '-', '.', , '/', '\\', '@', '[', ']', '(', ')', '#', '%', '=']);

  static SmartsRawAlphabet = new Set([
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '!', '#', '$', '&', '(', ')', '*', '+', ',', '-', '.', ':', ';', '=', '@', '~', '[', ']',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
    'N', 'O', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
    'n', 'o', 'p', 'r', 's', 't', 'u', 'v', 'y',
  ]);

  /** @param s {String} - string to check
   * @returns {boolean} */
  static isHelm(s) {
    return s.startsWith('PEPTIDE1{') || s.startsWith('CHEM1{') || s.startsWith('BLOB1{') ||
      s.startsWith('RNA1{') || s.startsWith('DNA1{');
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMacromolecule(col) {
    // To collect alphabet freq three strategies can be used:
    // as chars, as fasta (single or within square brackets), as with the separator.
    if (
      !(col.categories.length == 1 && !col.categories[0]) && // TODO: Remove with tests for single empty category value
      DG.Detector.sampleCategories(col, (s) => BioPackageDetectors.isHelm(s), 1)
    ) {
      col.setTag(DG.TAGS.UNITS, 'HELM');
      return DG.SEMTYPE.MACROMOLECULE;
    }

    const decoyAlphabets = [
      ['SMILES', BioPackageDetectors.SmilesRawAlphabet, 0.30],
      ['SMARTS', BioPackageDetectors.SmartsRawAlphabet, 0.45],
    ];

    const candidateAlphabets = [
      ['PT', BioPackageDetectors.PeptideFastaAlphabet, 0.55],
      ['DNA', BioPackageDetectors.DnaFastaAlphabet, 0.55],
      ['RNA', BioPackageDetectors.RnaFastaAlphabet, 0.55],
    ];

    // Check for url column, maybe it is too heavy check
    const isUrlCheck = (s) => {
      let res = true;
      try {
        const url = new URL(s);
        res = true;
      } catch {
        res = false;
      }
      return res;
    };
    const isUrl = DG.Detector.sampleCategories(col, isUrlCheck, 1);
    if (isUrl) return null;

    // TODO: Detect HELM sequence
    // TODO: Lazy calculations could be helpful for performance and convenient for expressing classification logic.
    const statsAsChars = BioPackageDetectors.getStats(col, 5, BioPackageDetectors.splitterAsChars);
    // if (Object.keys(statsAsChars.freq).length === 0) return;

    const decoy = BioPackageDetectors.detectAlphabet(statsAsChars.freq, decoyAlphabets, null);
    if (decoy != 'UN') return null;

    if (statsAsChars.sameLength) {
      if (Object.keys(statsAsChars.freq).length > 0) { // require non empty alphabet
        const alphabet = BioPackageDetectors.detectAlphabet(statsAsChars.freq, candidateAlphabets, '-');
        if (alphabet === 'UN') return null;

        //const units = `fasta:SEQ.MSA:${alphabet}`;
        const units = 'fasta';
        col.setTag(DG.TAGS.UNITS, units);
        col.setTag('aligned', 'SEQ.MSA');
        col.setTag('alphabet', alphabet);
        return DG.SEMTYPE.MACROMOLECULE;
      }
    } else {
      const separator = BioPackageDetectors.detectSeparator(statsAsChars.freq);
      const gapSymbol = separator ? '' : '-';
      const splitter = separator ? BioPackageDetectors.getSplitterWithSeparator(separator) : BioPackageDetectors.splitterAsFasta;

      const stats = BioPackageDetectors.getStats(col, 5, splitter);
      // Empty monomer alphabet is not allowed
      if (Object.keys(stats.freq).length === 0) return null;
      // Long monomer names for sequences with separators have constraints
      if (separator && BioPackageDetectors.checkForbiddenWithSeparators(stats.freq)) return null;

      const format = separator ? 'separator' : 'fasta';
      const seqType = stats.sameLength ? 'SEQ.MSA' : 'SEQ';

      // TODO: If separator detected, then extra efforts to detect alphabet are allowed.
      const alphabet = BioPackageDetectors.detectAlphabet(stats.freq, candidateAlphabets, gapSymbol);

      // const forbidden = BioPackageDetectors.checkForbiddenWoSeparator(stats.freq);
      if (separator || alphabet != 'UN') {
        //const units = `${format}:${seqType}:${alphabet}`;
        col.setTag(DG.TAGS.UNITS, format);
        col.setTag('aligned', seqType);
        col.setTag('alphabet', alphabet);
        if (separator) col.setTag('separator', separator);
        return DG.SEMTYPE.MACROMOLECULE;
      }
    }
  }

  /** Detects the most frequent char with a rate of at least 0.15 of others in sum.
   * Does not use any splitting strategies, estimates just by single characters.
   * */
  static detectSeparator(freq) {
    // To detect a separator we analyse col's sequences character frequencies.
    // If there is an exceptionally frequent symbol, then we will call it the separator.
    // The most frequent symbol should occur with a rate of at least 0.15
    // of all other symbols in sum to be called the separator.

    // !!! But there is a caveat because exceptionally frequent char can be a gap symbol in MSA.
    // !!! What is the difference between the gap symbol and separator symbol in stats terms?
    // const noSeparatorRe = /[a-z\d]+$/i;
    const noSeparatorChemRe = /[HBCNOFPSKVYI]/i; // Mendeleev's periodic table single char elements
    const noSeparatorAlphaDigitRe = /[\dA-Z,& _]/i; // ..., comma, ampersand, space, underscore
    const noSeparatorBracketsRe = /[\[\]()<>{}]/i;
    const cleanFreq = Object.assign({}, ...Object.entries(freq)
      .filter(([m, f]) =>
        !noSeparatorChemRe.test(m) && !noSeparatorAlphaDigitRe.test(m) && !noSeparatorBracketsRe.test(m) &&
        !BioPackageDetectors.PeptideFastaAlphabet.has(m) &&
        !BioPackageDetectors.DnaFastaAlphabet.has(m))
      .map(([m, f]) => ({[m]: f})));
    if (Object.keys(cleanFreq).length == 0) return null;

    const maxFreq = Math.max(...Object.values(cleanFreq));

    const sep = Object.entries(freq).find(([k, v]) => v === maxFreq)[0];
    const sepFreq = freq[sep];
    const otherSumFreq = Object.entries(freq).filter((kv) => kv[0] !== sep)
      .map((kv) => kv[1]).reduce((pSum, a) => pSum + a, 0);
    const freqThreshold = 3.5 * (1 / Object.keys(freq).length);
    return sepFreq / otherSumFreq > freqThreshold ? sep : null;
  }

  /** With a separator, spaces are nor allowed in monomer names.
   * The monomer name/label cannot contain digits only.
   */
  static checkForbiddenWithSeparators(freq) {
    const forbiddenRe = /[ ]|^\d+$/i;
    return Object.keys(freq).filter((m) => forbiddenRe.test(m)).length > 0;
  }

  // /** Without a separator, special symbols or digits are not allowed as monomers. */
  // static checkForbiddenWoSeparator(freq) {
  //   const forbiddenRe = /[\d!@#$%^&*()_+\-=\[\]{};':"\\|,.<>\/?]/i;
  //   return Object.keys(freq).filter((m) => forbiddenRe.test(m)).length > 0;
  // }

  /** Stats of sequences with specified splitter func, returns { freq, sameLength } */
  static getStats(seqCol, minLength, splitter) {
    const freq = {};
    let sameLength = true;
    let firstLength = null;

    for (const seq of seqCol.categories) {
      const mSeq = splitter(seq);

      if (firstLength == null) {
        firstLength = mSeq.length;
      } else if (mSeq.length !== firstLength) {
        sameLength = false;
      }

      if (mSeq.length > minLength) {
        for (const m of mSeq) {
          if (!(m in freq)) {
            freq[m] = 0;
          }
          freq[m] += 1;
        }
      }
    }
    return {freq: freq, sameLength: sameLength};
  }

  /** Detects alphabet for freq by freq similarity to alphabet monomer set.
   * @param freq       frequencies of monomers in sequence set
   * @param candidates  an array of pairs [name, monomer set]
   * */
  static detectAlphabet(freq, candidates, gapSymbol) {
    const candidatesSims = candidates.map((c) => {
      const sim = BioPackageDetectors.getAlphabetSimilarity(freq, c[1], gapSymbol);
      return [c[0], c[1], c[2], freq, sim];
    });

    let alphabetName;
    const maxSim = Math.max(...candidatesSims.map((cs) => cs[4] > cs[2] ? cs[4] : -1));
    if (maxSim > 0) {
      const sim = candidatesSims.find((cs) => cs[4] == maxSim);
      alphabetName = sim[0];
    } else {
      alphabetName = 'UN';
    }
    return alphabetName;
  }

  static getAlphabetSimilarity(freq, alphabet, gapSymbol) {
    const keys = new Set([...new Set(Object.keys(freq)), ...alphabet]);
    keys.delete(gapSymbol);

    const freqA = [];
    const alphabetA = [];
    for (const m of keys) {
      freqA.push(m in freq ? freq[m] : 0);
      alphabetA.push(alphabet.has(m) ? 10 : -20 /* penalty for character outside alphabet set*/);
    }
    /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
    const cos = BioPackageDetectors.vectorDotProduct(freqA, alphabetA) / (BioPackageDetectors.vectorLength(freqA) * BioPackageDetectors.vectorLength(alphabetA));
    return cos;
  }

  static vectorLength(v) {
    let sqrSum = 0;
    for (let i = 0; i < v.length; i++) {
      sqrSum += v[i] * v[i];
    }
    return Math.sqrt(sqrSum);
  }

  static vectorDotProduct(v1, v2) {
    if (v1.length != v2.length) {
      throw Error('The dimensionality of the vectors must match');
    }
    let prod = 0;
    for (let i = 0; i < v1.length; i++) {
      prod += v1[i] * v2[i];
    }
    return prod;
  }

  /** For trivial checks split by single chars*/
  static splitterAsChars(seq) {
    return seq.split('');
  }

  static getSplitterWithSeparator(separator) {
    return function(seq) {
      return seq.split(separator);
    };
  }

  // Multichar monomer names in square brackets, single char monomers or gap symbol
  static monomerRe = /\[(\w+)\]|(\w)|(-)/g;

  /** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol. */
  static splitterAsFasta(seq) {
    const res = wu(seq.toString().matchAll(BioPackageDetectors.monomerRe)).map((ma) => {
      let mRes;
      const m = ma[0];
      if (m.length > 1) {
        if (m in BioPackageDetectors.aaSynonyms) {
          mRes = BioPackageDetectors.aaSynonyms[m];
        } else {
          mRes = '';
          console.debug(`Long monomer '${m}' has not a short synonym.`);
        }
      } else {
        mRes = m;
      }
      return mRes;
    }).toArray();

    return res;
  }

  /** Only some of the synonyms. These were obtained from the clustered oligopeptide dataset. */
  static aaSynonyms = {
    '[MeNle]': 'L', // Nle - norleucine
    '[MeA]': 'A', '[MeG]': 'G', '[MeF]': 'F',
  };
}
