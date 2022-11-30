/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 *
 * TODO: Use detectors from WebLogo pickUp.. methods
 */

const SEQ_SAMPLE_LIMIT = 100;
const SEQ_SAMPLE_LENGTH_LIMIT = 500;

/** enum type to simplify setting "user-friendly" notation if necessary */
const NOTATION = {
  FASTA: 'fasta',
  SEPARATOR: 'separator',
  HELM: 'helm',
};

const ALPHABET = {
  DNA: 'DNA',
  RNA: 'RNA',
  PT: 'PT',
  UN: 'UN',
};

const ALIGNMENT = {
  SEQ_MSA: 'SEQ.MSA',
  SEQ: 'SEQ',
};

/** Class for handling notation units in Macromolecule columns */
const UnitsHandler = {
  TAGS: {
    aligned: 'aligned',
    alphabet: 'alphabet',
    alphabetSize: '.alphabetSize',
    alphabetIsMultichar: '.alphabetIsMultichar',
    separator: 'separator',
  },
};

const isUrlRe = /[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}\b([-a-zA-Z0-9()@:%_\+.~#?&//=]*)?/i;

class BioPackageDetectors extends DG.Package {

  PeptideFastaAlphabet = new Set([
    'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
    'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
    'MeNle', 'MeA', 'MeG', 'MeF',
  ]);

  DnaFastaAlphabet = new Set(['A', 'C', 'G', 'T']);

  RnaFastaAlphabet = new Set(['A', 'C', 'G', 'U']);

  SmilesRawAlphabet = new Set([
    'A', 'B', 'C', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'R', 'S', 'Z',
    'a', 'c', 'e', 'g', 'i', 'l', 'n', 'o', 'r', 's', 't', 'u',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '+', '-', '.', , '/', '\\', '@', '[', ']', '(', ')', '#', '%', '=']);

  SmartsRawAlphabet = new Set([
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '!', '#', '$', '&', '(', ')', '*', '+', ',', '-', '.', ':', ';', '=', '@', '~', '[', ']',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
    'N', 'O', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
    'n', 'o', 'p', 'r', 's', 't', 'u', 'v', 'y',
  ]);

  /** @param s {String} - string to check
   * @returns {boolean} */
  isHelm(s) {
    return s.startsWith('PEPTIDE1{') || s.startsWith('CHEM1{') || s.startsWith('BLOB1{') ||
      s.startsWith('RNA1{') || s.startsWith('DNA1{');
  }

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMacromolecule(col) {
    const t1 = Date.now();
    try {
      // Fail early
      if (col.type !== DG.TYPE.STRING) return null;

      const categoriesSample = col.categories.length < SEQ_SAMPLE_LIMIT ? col.categories :
        this.sample(col.categories, SEQ_SAMPLE_LIMIT);

      // To collect alphabet freq three strategies can be used:
      // as chars, as fasta (single or within square brackets), as with the separator.
      if (
        !(col.categories.length == 1 && !col.categories[0]) && // TODO: Remove with tests for single empty category value
        DG.Detector.sampleCategories(col, (s) => this.isHelm(s), 1, SEQ_SAMPLE_LIMIT)
      ) {
        const statsAsHelm = this.getStats(categoriesSample, 2,
          this.getSplitterAsHelm(SEQ_SAMPLE_LENGTH_LIMIT));
        col.setTag(DG.TAGS.UNITS, NOTATION.HELM);

        // alphabetSize calculated on (sub)sample of data is incorrect
        // const alphabetSize = Object.keys(statsAsHelm.freq).length;
        const alphabetIsMultichar = Object.keys(statsAsHelm.freq).some((m) => m.length > 1);
        // col.setTag(UnitsHandler.TAGS.alphabetSize, alphabetSize.toString());
        col.setTag(UnitsHandler.TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');

        return DG.SEMTYPE.MACROMOLECULE;
      }

      const decoyAlphabets = [
        ['SMILES', this.SmilesRawAlphabet, 0.30],
        ['SMARTS', this.SmartsRawAlphabet, 0.43],
      ];

      const candidateAlphabets = [
        [ALPHABET.PT, this.PeptideFastaAlphabet, 0.50],
        [ALPHABET.DNA, this.DnaFastaAlphabet, 0.55],
        [ALPHABET.RNA, this.RnaFastaAlphabet, 0.55],
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
        // return isUrlRe.test(s);
      };
      const isUrl = categoriesSample.every((v) => { return !v || isUrlCheck(v); });
      if (isUrl) return null;

      // TODO: Detect HELM sequence
      // TODO: Lazy calculations could be helpful for performance and convenient for expressing classification logic.
      const statsAsChars = this.getStats(categoriesSample, 5,
        this.getSplitterAsChars(SEQ_SAMPLE_LENGTH_LIMIT));
      // Empty statsAsShars.freq alphabet means no strings of enough length presented in the data
      if (Object.keys(statsAsChars.freq).length === 0) return null;

      const decoy = this.detectAlphabet(statsAsChars.freq, decoyAlphabets, null);
      if (decoy != ALPHABET.UN) return null;

      const separator = this.detectSeparator(statsAsChars.freq);
      const units = separator ? NOTATION.SEPARATOR : NOTATION.FASTA;
      const gapSymbol = separator ? '' : '-';
      const splitter = separator ? this.getSplitterWithSeparator(separator, SEQ_SAMPLE_LENGTH_LIMIT) :
        this.getSplitterAsFasta(SEQ_SAMPLE_LENGTH_LIMIT);

      if (statsAsChars.sameLength) {
        const stats = this.getStats(categoriesSample, 5, splitter);
        const alphabet = this.detectAlphabet(stats.freq, candidateAlphabets, '-');
        if (alphabet === ALPHABET.UN) return null;

        col.setTag(DG.TAGS.UNITS, units);
        if (separator) col.setTag(UnitsHandler.TAGS.separator, separator);
        col.setTag(UnitsHandler.TAGS.aligned, ALIGNMENT.SEQ_MSA);
        col.setTag(UnitsHandler.TAGS.alphabet, alphabet);
        return DG.SEMTYPE.MACROMOLECULE;
      } else {
        const stats = this.getStats(categoriesSample, 5, splitter);
        // Empty monomer alphabet is not allowed
        if (Object.keys(stats.freq).length === 0) return null;
        // Long monomer names for sequences with separators have constraints
        if (separator && this.checkForbiddenWithSeparators(stats.freq)) return null;

        const aligned = stats.sameLength ? ALIGNMENT.SEQ_MSA : ALIGNMENT.SEQ;

        // TODO: If separator detected, then extra efforts to detect alphabet are allowed.
        const alphabet = this.detectAlphabet(stats.freq, candidateAlphabets, gapSymbol);

        // const forbidden = this.checkForbiddenWoSeparator(stats.freq);
        if (separator || alphabet != 'UN') {
          col.setTag(DG.TAGS.UNITS, units);
          if (separator) col.setTag(UnitsHandler.TAGS.separator, separator);
          col.setTag(UnitsHandler.TAGS.aligned, aligned);
          col.setTag(UnitsHandler.TAGS.alphabet, alphabet);
          if (alphabet === ALPHABET.UN) {
            // alphabetSize calculated on (sub)sample of data is incorrect
            // const alphabetSize = Object.keys(stats.freq).length;
            const alphabetIsMultichar = Object.keys(stats.freq).some((m) => m.length > 1);
            // col.setTag(UnitsHandler.TAGS.alphabetSize, alphabetSize.toString());
            col.setTag(UnitsHandler.TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');
          }
          return DG.SEMTYPE.MACROMOLECULE;
        }
      }
    } finally {
      const t2 = Date.now();
      console.debug('Bio: detectMacromolecule() ' + `ET = ${t2 - t1} ms.`);
    }
  }

  /** Detects the most frequent char with a rate of at least 0.15 of others in sum.
   * Does not use any splitting strategies, estimates just by single characters.
   * */
  detectSeparator(freq) {
    // To detect a separator we analyse col's sequences character frequencies.
    // If there is an exceptionally frequent symbol, then we will call it the separator.
    // The most frequent symbol should occur with a rate of at least 0.15
    // of all other symbols in sum to be called the separator.

    // !!! But there is a caveat because exceptionally frequent char can be a gap symbol in MSA.
    // !!! What is the difference between the gap symbol and separator symbol in stats terms?
    // const noSeparatorRe = /[a-z\d]+$/i;
    const noSeparatorChemRe = /[HBCNOFPSKVYI]/i; // Mendeleev's periodic table single char elements
    const noSeparatorAlphaDigitRe = /[\dA-Z,& _\r\n]/i; // ..., comma, ampersand, space, underscore, CR, LF
    const noSeparatorBracketsRe = /[\[\]()<>{}]/i;
    const cleanFreq = Object.assign({}, ...Object.entries(freq)
      .filter(([m, f]) =>
        !noSeparatorChemRe.test(m) && !noSeparatorAlphaDigitRe.test(m) && !noSeparatorBracketsRe.test(m) &&
        !this.PeptideFastaAlphabet.has(m) &&
        !this.DnaFastaAlphabet.has(m))
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
  checkForbiddenWithSeparators(freq) {
    const forbiddenRe = /[ ]|^\d+$/i;
    return Object.keys(freq).filter((m) => forbiddenRe.test(m)).length > 0;
  }

  // /** Without a separator, special symbols or digits are not allowed as monomers. */
  //  checkForbiddenWoSeparator(freq) {
  //   const forbiddenRe = /[\d!@#$%^&*()_+\-=\[\]{};':"\\|,.<>\/?]/i;
  //   return Object.keys(freq).filter((m) => forbiddenRe.test(m)).length > 0;
  // }

  /** Stats of sequences with specified splitter func, returns { freq, sameLength } */
  getStats(values, minLength, splitter) {
    const freq = {};
    let sameLength = true;
    let firstLength = null;

    for (const seq of values) {
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
  detectAlphabet(freq, candidates, gapSymbol) {
    const candidatesSims = candidates.map((c) => {
      const sim = this.getAlphabetSimilarity(freq, c[1], gapSymbol);
      return [c[0], c[1], c[2], freq, sim];
    });

    let alphabetName;
    const maxSim = Math.max(...candidatesSims.map((cs) => cs[4] > cs[2] ? cs[4] : -1));
    if (maxSim > 0) {
      const sim = candidatesSims.find((cs) => cs[4] == maxSim);
      alphabetName = sim[0];
    } else {
      alphabetName = ALPHABET.UN;
    }
    return alphabetName;
  }

  getAlphabetSimilarity(freq, alphabet, gapSymbol) {
    const keys = new Set([...new Set(Object.keys(freq)), ...alphabet]);
    keys.delete(gapSymbol);

    const freqA = [];
    const alphabetA = [];
    for (const m of keys) {
      freqA.push(m in freq ? freq[m] : 0);
      alphabetA.push(alphabet.has(m) ? 10 : -20 /* penalty for character outside alphabet set*/);
    }
    /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
    const cos = this.vectorDotProduct(freqA, alphabetA) / (this.vectorLength(freqA) * this.vectorLength(alphabetA));
    return cos;
  }

  vectorLength(v) {
    let sqrSum = 0;
    for (let i = 0; i < v.length; i++) {
      sqrSum += v[i] * v[i];
    }
    return Math.sqrt(sqrSum);
  }

  vectorDotProduct(v1, v2) {
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
  getSplitterAsChars(lengthLimit) {
    return function(seq) {
      return seq.split('', lengthLimit);
    }.bind(this);
  }

  getSplitterWithSeparator(separator, lengthLimit) {
    return function(seq) {
      // if (!!lengthLimit) {
      //   const res = new Array(lengthLimit);
      //   let pos = 0, count = 0;
      //   while (pos < seq.length && count < lengthLimit) {
      //     const newPos = seq.indexOf(separator, pos);
      //     res[count] = seq.substring(pos, newPos);
      //     count++;
      //     pos = newPos;
      //   }
      //
      //   return res.slice(0, count);
      // } else {
      return seq.split(separator, lengthLimit);
      // }
    }.bind(this);
  }

  // Multichar monomer names in square brackets, single char monomers or gap symbol
  monomerRe = /\[(\w+)\]|(\w)|(-)/g;

  /** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol. */
  getSplitterAsFasta(lengthLimit) {
    return function(seq) {
      const res = wu(seq.toString().matchAll(this.monomerRe))
        .take(lengthLimit)
        .map((ma) => {
          let mRes;
          const m = ma[0];
          if (m.length > 1) {
            mRes = ma[1];
          } else {
            mRes = m;
          }
          return mRes;
        }).toArray();

      return res;
    }.bind(this);
  }

  /** Only some of the synonyms. These were obtained from the clustered oligopeptide dataset. */
  aaSynonyms = {
    '[MeNle]': 'L', // Nle - norleucine
    '[MeA]': 'A', '[MeG]': 'G', '[MeF]': 'F',
  };

  helmRe = /(PEPTIDE1|DNA1|RNA1)\{([^}]+)}/g;
  helmPp1Re = /\[([^\[\]]+)]/g;

  /** Splits Helm string to monomers, but does not replace monomer names to other notation (e.g. for RNA). */
  getSplitterAsHelm(lengthLimit) {
    return function(seq) {
      this.helmRe.lastIndex = 0;
      const ea = this.helmRe.exec(seq.toString());
      const inSeq = ea ? ea[2] : null;

      const mmPostProcess = (mm) => {
        this.helmPp1Re.lastIndex = 0;
        const pp1M = this.helmPp1Re.exec(mm);
        if (pp1M && pp1M.length >= 2) {
          return pp1M[1];
        } else {
          return mm;
        }
      };

      const mmList = inSeq ? inSeq.split('.') : [];
      const mmListRes = mmList.map(mmPostProcess);
      return mmListRes;
    }.bind(this);
  }

  sample(src, n) {
    if (src.length < n) {
      throw new Error('Sample source is less than n requested.');
    }

    const idxSet = new Set();
    while (idxSet.size < n) {
      const idx = Math.floor(Math.random() * src.length);
      if (!idxSet.has(idx)) {
        idxSet.add(idx);
      }
    }

    return [...idxSet].map((idx) => src[idx]);
  }
}
