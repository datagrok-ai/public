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

  static semType = 'MACROMOLECULE';

  static Units = {
    FastaSeqPt: 'fasta:SEQ:PT', FastaSeqNt: 'fasta:SEQ:NT', FastaMsaPt: 'fasta:MSA:PT', FastaMsaNt: 'fasta:MSA:NT',
  };

  static AminoacidsFastaAlphabet = new Set([
    'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
    'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
  ]);

  static NucleotidesFastaAlphabet = new Set(['A', 'C', 'G', 'T']);

  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMacromolecule(col) {
    // To collect alphabet freq three strategies can be used:
    // as chars, as fasta (single or within square brackets), as with the separator.

    const alphabetCandidates = [
      ['NT', BioPackageDetectors.NucleotidesFastaAlphabet],
      ['PT', BioPackageDetectors.AminoacidsFastaAlphabet],
    ];

    // TODO: Detect HELM sequence
    // TODO: Lazy calculations could be helpful for performance and convenient for expressing classification logic.
    const statsAsChars = BioPackageDetectors.getStats(col, 5, BioPackageDetectors.splitterAsChars);
    if (statsAsChars.sameLength) {
      const alphabet = BioPackageDetectors.detectAlphabet(statsAsChars.freq, alphabetCandidates, '-');
      const units = `fasta:SEQ.MSA:${alphabet}`;
      col.setTag(DG.TAGS.UNITS, units);
      return BioPackageDetectors.semType;
    } else {
      const sep = BioPackageDetectors.detectSeparator(statsAsChars.freq);
      const gapSymbol = sep ? '' : '-';
      const splitter = sep ? BioPackageDetectors.getSplitterWithSeparator(sep) : BioPackageDetectors.splitterAsFasta;
      const stats = BioPackageDetectors.getStats(col, 5, splitter);

      const format = sep ? 'separator' : 'fasta';
      const seqType = stats.sameLength ? 'SEQ.MSA' : 'SEQ';

      // TODO: If separator detected, then extra efforts to detect alphabet are allowed.
      const alphabet = BioPackageDetectors.detectAlphabet(stats.freq, alphabetCandidates, gapSymbol);

      const units = `${format}:${seqType}:${alphabet}`;
      col.setTag(DG.TAGS.UNITS, units);
      return BioPackageDetectors.semType;
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

    const maxFreq = Math.max(...Object.values(freq));
    const sep = Object.entries(freq).find((kv) => kv[1] == maxFreq)[0];
    const sepFreq = freq[sep];
    const otherSumFreq = Object.entries(freq).filter((kv) => kv[0] !== sep)
      .map((kv) => kv[1]).reduce((pSum, a) => pSum + a, 0);
    const freqThreshold = 3.5 * (1 / Object.keys(freq).length);
    return sepFreq / otherSumFreq > freqThreshold ? sep : null;
  }

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
      return [c[0], c[1], freq, sim];
    });

    let alphabetName;
    const maxSim = Math.max(...candidatesSims.map((cs) => cs[3]));
    if (maxSim > 0.65) {
      const sim = candidatesSims.find((cs) => cs[3] == maxSim);
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
      alphabetA.push(alphabet.has(m) ? 1 : 0);
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

  static getSplitterWithSeparator(sep) {
    return function(seq) {
      return seq.split(sep);
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
