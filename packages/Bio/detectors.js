/* eslint-disable max-lines-per-function */
/* eslint-disable max-lines */
'use strict';
/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 *
 * TODO: Use detectors from WebLogo pickUp.. methods
 */
// eslint-disable-next-line max-lines

const SEQ_SAMPLE_LIMIT = 100;
const SEQ_SAMPLE_LENGTH_LIMIT = 100;

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

const SeqTemps = {
  seqHandler: `seq-handler`,
  notationProvider: `seq-handler.notation-provider`,
};


/** Class for handling notation units in Macromolecule columns */
const SeqHandler = {
  TAGS: {
    aligned: 'aligned',
    alphabet: 'alphabet',
    alphabetSize: '.alphabetSize',
    alphabetIsMultichar: '.alphabetIsMultichar',
    separator: 'separator',
  },
};

const isUrlRe = /[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}\b([-a-zA-Z0-9()@:%_\+.~#?&//=]*)?/i;

class LoggerWrapper {
  constructor(_package, logger, componentName) {
    this.package = _package;
    this.logger = logger;
    this.componentName = componentName;
    this.debugEnabled = false;
  }

  debug(message, params) {
    if (!this.debugEnabled) return;
    this.logger.debug(message, params);
  }

  error(message, params, stackTrace) {
    this.logger.error(message, params, stackTrace);
  }
}

class BioPackageDetectors extends DG.Package {
  static objCounter = -1;
  objId = ++BioPackageDetectors.objCounter;

  constructor() {
    super();

    this.forbiddenMulticharAll = ' .:';
    this.forbiddenMulticharFirst = ']' + this.forbiddenMulticharAll;
    this.forbiddenMulticharMiddle = '][' + this.forbiddenMulticharAll;
    this.forbiddenMulticharLast = '[' + this.forbiddenMulticharAll;

    // replace super._logger
    this._logger = new LoggerWrapper(this, this.logger, 'detectors');
  }

  /** Parts of the column name required in the column's name under the detector. It must be in lowercase. */
  likelyColNamePartList = ['seq', 'msa', 'dna', 'rna', 'fasta', 'helm', 'sense', 'protein'];

  peptideFastaAlphabet = new Set([
    'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
    'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
    'MeNle', 'MeA', 'MeG', 'MeF',
  ]);

  dnaFastaAlphabet = new Set(['A', 'C', 'G', 'T']);

  rnaFastaAlphabet = new Set(['A', 'C', 'G', 'U']);

  numbersRawAlphabet = new Set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']);

  smilesRawAlphabet = new Set([
    'C', 'F', 'H', 'N', 'O', 'P', 'S', 'B', /**/'A', 'E', 'I', 'K', 'L', 'M', 'R', 'Z',
    'c', 'n', 'o', 's', /**/'a', 'e', 'g', 'i', 'l', 'r', 't', 'u',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    '+', '-', '.', , '/', '\\', '@', '[', ']', '(', ')', '#', '%', '=']);

  smartsRawAlphabet = new Set([
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

  //name: detectMacromoleculeEnableStore
  //output: object result
  detectMacromoleculeEnableStore() {
    return window.$detectMacromoleculeStore = {last: null};
  }

  /** Returns last object (stores it if enabled earlier). */
  detectMacromoleculeStoreLast() {
    const last = {};
    if (window.$detectMacromoleculeStore) window.$detectMacromoleculeStore.last = last;
    return last;
  }

  /** Detector MUST NOT be async, causes error:
   *  Concurrent modification during iteration: Instance of 'JSArray<Column>'.
   */
  //tags: semTypeDetector
  //input: column col
  //output: string semType
  detectMacromolecule(col) {
    const tableName = col.dataFrame ? col.dataFrame.name : null;
    this.logger.debug(`Bio: detectMacromolecule( table: ${tableName}.${col.name} ), start`);
    const t1 = Date.now();
    try {
      const last = this.detectMacromoleculeStoreLast();
      const colName = col.name;
      const colNameLikely = this.likelyColNamePartList.some(
        (requiredColNamePart) => colName.toLowerCase().includes(requiredColNamePart));
      const seqMinLength = colNameLikely ? 7 : 10;
      const maxBadRatio = colNameLikely ? 0.05 : 0.005;

      // Fail early
      if (col.type !== DG.TYPE.STRING) {
        last.rejectReason = `The column must be of type '${DG.TYPE.STRING}'.`;
        return null;
      }

      const categoriesSample = [...new Set((col.length < SEQ_SAMPLE_LIMIT ?
        wu.count(0).take(Math.min(SEQ_SAMPLE_LIMIT, col.length)).map((rowI) => col.get(rowI)) :
        this.sample(col, SEQ_SAMPLE_LIMIT))
        .map((seq) => !!seq ? seq.substring(0, SEQ_SAMPLE_LENGTH_LIMIT * 5) : '')
        .filter((seq) => seq.length !== 0/* skip empty values for detector */),
      )];
      last.categoriesSample = categoriesSample;

      // To collect alphabet freq three strategies can be used:
      // as chars, as fasta (single or within square brackets), as with the separator.
      if (
        !(col.categories.length === 1 && !col.categories[0]) && // TODO: Remove with tests for single empty category
        DG.Detector.sampleCategories(col, (s) => this.isHelm(s), 1, SEQ_SAMPLE_LIMIT)
      ) {
        const statsAsHelm = this.getStats(categoriesSample, 2,
          this.getSplitterAsHelm(SEQ_SAMPLE_LENGTH_LIMIT));
        col.meta.units = NOTATION.HELM;

        // alphabetSize calculated on (sub)sample of data is incorrect
        // const alphabetSize = Object.keys(statsAsHelm.freq).length;
        const alphabetIsMultichar = Object.keys(statsAsHelm.freq).some((m) => m.length > 1);
        // col.setTag(SeqHandler.TAGS.alphabetSize, alphabetSize.toString());
        col.setTag(SeqHandler.TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');

        col.setTag(DG.TAGS.CELL_RENDERER, 'helm');
        return DG.SEMTYPE.MACROMOLECULE;
      }

      const decoyAlphabets = [
        ['NUMBERS', this.numbersRawAlphabet, 0.25, undefined],
        ['SMILES', this.smilesRawAlphabet, 0.25, (seq) => seq.replaceAll()],
        ['SMARTS', this.smartsRawAlphabet, 0.45, undefined],
      ];

      const candidateAlphabets = [
        [ALPHABET.PT, this.peptideFastaAlphabet, 0.50],
        [ALPHABET.DNA, this.dnaFastaAlphabet, 0.55],
        [ALPHABET.RNA, this.rnaFastaAlphabet, 0.55],
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
      const isUrl = categoriesSample.every((v) => !v || isUrlCheck(v));
      if (isUrl) {
        last.rejectReason = 'URL detected.';
        return null;
      }

      // TODO: Detect HELM sequence
      // TODO: Lazy calculations could be helpful for performance and convenient for expressing classification logic.
      const statsAsChars = this.getStats(categoriesSample, seqMinLength,
        this.getSplitterAsChars(SEQ_SAMPLE_LENGTH_LIMIT));
      // Empty statsAsShars.freq alphabet means no strings of enough length presented in the data
      if (Object.keys(statsAsChars.freq).length === 0) {
        last.rejectReason = 'Monomer set (alphabet) is empty.';
        return null;
      }

      const decoy = this.detectAlphabet(statsAsChars.freq, decoyAlphabets, null, colNameLikely ? -0.05 : 0);
      if (decoy !== ALPHABET.UN) {
        last.rejectReason = `Decoy alphabet '${decoy}' detected.`;
        return null;
      }

      const separator = this.detectSeparator(statsAsChars.freq, categoriesSample, seqMinLength);
      const checkForbiddenSeparatorRes = this.checkForbiddenSeparator(separator);
      if (checkForbiddenSeparatorRes) {
        last.rejectReason = `Separator '${separator}' is forbidden.`;
        return null;
      }

      const units = separator ? NOTATION.SEPARATOR : NOTATION.FASTA;
      const gapSymbol = separator ? '' : '-';
      const splitter = separator ? this.getSplitterWithSeparator(separator, SEQ_SAMPLE_LENGTH_LIMIT) :
        this.getSplitterAsFasta(SEQ_SAMPLE_LENGTH_LIMIT);

      if (statsAsChars.sameLength && !separator &&
        !(['[', ']'].some((c) => c in statsAsChars.freq)) // not fasta ext notation
      ) { // MSA FASTA single character
        const stats = this.getStats(categoriesSample, seqMinLength, splitter);
        const alphabet = this.detectAlphabet(stats.freq, candidateAlphabets, '-', colNameLikely ? 0.20 : 0);
        if (alphabet === ALPHABET.UN) {
          last.rejectReason = `MSA FASTA single character alphabet is not allowed to be 'UN'.`;
          return null;
        }

        col.meta.units = units;
        if (separator) col.setTag(SeqHandler.TAGS.separator, separator);
        col.setTag(SeqHandler.TAGS.aligned, ALIGNMENT.SEQ_MSA);
        col.setTag(SeqHandler.TAGS.alphabet, alphabet);
        if (alphabet === ALPHABET.UN) {
          const alphabetIsMultichar = Object.keys(stats.freq).some((m) => m.length > 1);
          col.setTag(SeqHandler.TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');
        }
        col.setTag(DG.TAGS.CELL_RENDERER, 'sequence');
        return DG.SEMTYPE.MACROMOLECULE;
      } else {
        // for fasta, we need to include every sequence
        const stats = this.getStats(categoriesSample, separator ? seqMinLength : 2, splitter);
        const alphabetIsMultichar = Object.keys(stats.freq).some((m) => m.length > 1);
        // Empty monomer alphabet is not allowed
        if (Object.keys(stats.freq).length === 0) {
          last.rejectReason = 'Monomer set (alphabet) is empty';
          return null;
        }
        // Single- and multi-char monomer names for sequences with separators have constraints
        if (units === NOTATION.SEPARATOR || (units === NOTATION.FASTA && alphabetIsMultichar)) {
          const badSymbol /*: string | null*/ = this.checkBadMultichar(stats.freq);
          if (badSymbol) {
            last.rejectReason = `Forbidden multi-char monomer: '${badSymbol}'.`;
            return null;
          }
        }
        const aligned = stats.sameLength ? ALIGNMENT.SEQ_MSA : ALIGNMENT.SEQ;

        // TODO: If separator detected, then extra efforts to detect alphabet are allowed.
        const alphabet = this.detectAlphabet(stats.freq, candidateAlphabets, gapSymbol, colNameLikely ? 0.15 : 0);
        if (units === NOTATION.FASTA && alphabet === ALPHABET.UN && !alphabetIsMultichar) {
          last.rejectReason = `FASTA single character alphabet is not allowed to be 'UN'.`;
          return null;
        }

        // const forbidden = this.checkForbiddenWoSeparator(stats.freq);
        col.meta.units = units;
        if (separator) col.setTag(SeqHandler.TAGS.separator, separator);
        col.setTag(SeqHandler.TAGS.aligned, aligned);
        col.setTag(SeqHandler.TAGS.alphabet, alphabet);
        if (alphabet === ALPHABET.UN) {
          // alphabetSize calculated on (sub)sample of data is incorrect
          const alphabetIsMultichar = Object.keys(stats.freq).some((m) => m.length > 1);
          col.setTag(SeqHandler.TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');
        }

        refineSeqSplitter(col, stats, separator).then(() => { });
        col.setTag(DG.TAGS.CELL_RENDERER, 'sequence');
        return DG.SEMTYPE.MACROMOLECULE;
      }
    } catch (err) {
      const errMsg = err instanceof Error ? err.message : err.toString();
      const errStack = err instanceof Error ? err.stack : undefined;
      const colTops = wu.count(0).take(Math.max(col.length, 4)).map((rowI) => col.get(rowI))
        .reduce((a, b) => a === undefined ? b : a + '\n' + b, undefined);
      this.logger.error(`Bio: detectMacromolecule( table: ${tableName}.${col.name} ), error:\n${errMsg}` +
        `${errStack ? '\n' + errStack : ''}` + `\n${colTops}`);
    } finally {
      // Prevent too much log spam
      // const t2 = Date.now();
      // this.logger.debug(`Bio: detectMacromolecule( table: ${tableName}.${col.name} ), ` + `ET = ${t2 - t1} ms.`);
    }
  }

  /** Detects the most frequent char with a rate of at least 0.15 of others in sum.
   * Does not use any splitting strategies, estimates just by single characters.
   * @param freq Dictionary of characters freqs
   * @param categoriesSample A string array of seqs sample
   * @param seqMinLength A threshold on min seq length for contributing to stats
   */
  detectSeparator(freq, categoriesSample, seqMinLength) {
    // To detect a separator we analyze col's sequences character frequencies.
    // If there is an exceptionally frequent symbol, then we will call it the separator.
    // The most frequent symbol should occur with a rate of at least 0.15
    // of all other symbols in sum to be called the separator.

    // !!! But there is a caveat because exceptionally frequent char can be a gap symbol in MSA.
    // !!! What is the difference between the gap symbol and separator symbol in stats terms?
    // const noSeparatorRe = /[a-z\d]+$/i;
    const noSeparatorChemRe = /[HBCNOFPSKVYI]/i; // Mendeleev's periodic table single char elements
    const noSeparatorAlphaDigitRe = /[\dA-Z]/i;
    const noSeparatorBracketsRe = /[\[\]()<>{}]/i;
    const cleanFreq = Object.assign({}, ...Object.entries(freq)
      .filter(([m, f]) =>
        !noSeparatorChemRe.test(m) && !noSeparatorAlphaDigitRe.test(m) && !noSeparatorBracketsRe.test(m) &&
        !this.peptideFastaAlphabet.has(m) &&
        !this.dnaFastaAlphabet.has(m))
      .map(([m, f]) => ({[m]: f})));
    if (Object.keys(cleanFreq).length === 0) return null;

    const maxFreq = Math.max(...Object.values(cleanFreq));

    const sep = Object.entries(freq).find(([k, v]) => v === maxFreq)[0];
    const sepFreq = freq[sep];
    const otherSumFreq = Object.entries(freq).filter((kv) => kv[0] !== sep)
      .map((kv) => kv[1]).reduce((pSum, a) => pSum + a, 0);

    // Splitter with separator test application
    const splitter = this.getSplitterWithSeparator(sep, SEQ_SAMPLE_LENGTH_LIMIT);
    const stats = this.getStats(categoriesSample, seqMinLength, splitter);
    const badSymbol = this.checkBadMultichar(stats.freq);
    if (badSymbol) return null;
    // TODO: Test for Gamma/Erlang distribution
    const totalMonomerCount = wu(Object.values(stats.freq)).reduce((sum, a) => sum + a, 0);
    const mLengthAvg = wu.entries(stats.freq)
      .reduce((sum, [m, c]) => sum + m.length * c, 0) / totalMonomerCount;
    const mLengthVarN = Math.sqrt(wu.entries(stats.freq)
      .reduce((sum, [m, c]) => sum + Math.pow(m.length - mLengthAvg, 2) * c, 0) / (totalMonomerCount - 1),
    ) / mLengthAvg;

    const sepRate = sepFreq / (sepFreq + otherSumFreq);
    const expSepRate = 1 / Object.keys(freq).length; // expected
    // const freqThreshold = (1 / (Math.log2(Object.keys(freq).length) + 2));

    return (sepRate / expSepRate > 2.2 && mLengthVarN < 0.8) ||
    (sepRate / expSepRate > 3.5) ? sep : null;
  }

  checkForbiddenSeparator(separator) {
    // comma, ampersand, space, underscore, CRLF, CR, LF
    // 2023-04-15: dot is allowed to allow Helm like separator in Helm MSA results (no Helm monomers contains dot)
    const forbiddenSepRe = /,|&| |_|\r\n|\r|\n/i;
    return forbiddenSepRe.test(separator);
  }

  /** Dots and colons are nor allowed in multichar monomer names (but space is allowed).
   * The monomer name/label cannot contain digits only (but single digit is allowed).
   */
  checkBadMultichar(freq) /* : string | null */ {
    for (const symbol of Object.keys(freq)) {
      if (symbol && !isNaN(symbol))
        return symbol; // performance evaluated better with RegExp

      const symbolLen = symbol.length;
      if (this.forbiddenMulticharFirst.includes(symbol[0]))
        return symbol;
      if (this.forbiddenMulticharLast.includes(symbol[symbolLen - 1]))
        return symbol;
      for (let cI = 1; cI < symbolLen - 1; ++cI) {
        const c = symbol[cI];
        if (this.forbiddenMulticharMiddle.includes(c))
          return symbol;
      }
      if (symbol.match(/^\d+\W+.*/))
        // symbols like '2,...' are forbidden
        // we require an alphabet character just after the leading digit(s)
        return symbol;
    }
    return null;
  }

  calcBad(freq, forbiddenSet) {
    let allCount = 0;
    let forbiddenCount = 0;
    for (const [m, count] of Object.entries(freq)) {
      if (forbiddenSet.has(m)) forbiddenCount += freq[m];
      allCount += freq[m];
    }
    return [forbiddenCount, allCount];
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
      const mSeq = !!seq ? splitter(seq) : [];

      if (firstLength === null)
        firstLength = mSeq.length;
      else if (mSeq.length !== firstLength)
        sameLength = false;

      if (mSeq.length >= minLength) {
        for (const m of mSeq) {
          if (!(m in freq)) freq[m] = 0;
          freq[m] += 1;
        }
      }
    }
    return {freq: freq, sameLength: sameLength};
  }

  /** Detects alphabet for freq by freq similarity to alphabet monomer set.
   * @param freq       frequencies of monomers in sequence set
   * @param candidates  an array of pairs [name, monomer set]
   * @param {boolean} colNameLikely The column name suggests the column is Macromolecule more likely
   */
  detectAlphabet(freq, candidates, gapSymbol, simAdj = 0) {
    const candidatesSims = candidates.map((c) => {
      const sim = this.getAlphabetSimilarity(freq, c[1], gapSymbol) + simAdj;
      return [c[0], c[1], c[2], freq, sim];
    });

    let alphabetName;
    const maxSim = Math.max(...candidatesSims.map((cs) => cs[4] > cs[2] ? cs[4] : -1));
    if (maxSim > 0) {
      const sim = candidatesSims.find((cs) => cs[4] === maxSim);
      alphabetName = sim[0];
    } else
      alphabetName = ALPHABET.UN;
    return alphabetName;
  }

  getAlphabetSimilarity(freq, alphabet, gapSymbol) {
    const keys = new Set([...new Set(Object.keys(freq)), ...alphabet]);
    keys.delete(gapSymbol);

    const freqSum = Object.values(freq).reduce((a, b) => a + b, 0);
    const freqA = [];
    const alphabetA = [];
    for (const m of keys) {
      freqA.push(m in freq ? freq[m] / freqSum : -0.001);
      alphabetA.push(alphabet.has(m) ? 10 : -20 /* penalty for character outside alphabet set*/);
    }
    /* There were a few ideas: chi-squared, pearson correlation (variance?), scalar product */
    const cos = this.vectorDotProduct(freqA, alphabetA) / (this.vectorLength(freqA) * this.vectorLength(alphabetA));
    return cos;
  }

  vectorLength(v) {
    let sqrSum = 0;
    for (let i = 0; i < v.length; i++)
      sqrSum += v[i] * v[i];
    return Math.sqrt(sqrSum);
  }

  vectorDotProduct(v1, v2) {
    if (v1.length !== v2.length)
      throw Error('The dimensionality of the vectors must match');

    let prod = 0;
    for (let i = 0; i < v1.length; i++)
      prod += v1[i] * v2[i];

    return prod;
  }

  /** For trivial checks split by single chars*/
  getSplitterAsChars(lengthLimit) {
    const resFunc = function(seq) {
      return seq.split('', lengthLimit);
    };
    resFunc.T = 'splitterAsChars';
    return resFunc;
  }

  getSplitterWithSeparator(separator, limit) {
    const resFunc = function(seq) {
      return !seq ? [] : seq.replaceAll('\"-\"', '').replaceAll('\'-\'', '').split(separator, limit);
    };
    resFunc.T = 'splitterWithSeparator';
    return resFunc;
  }

  // Multichar monomer names in square brackets, single char monomers or gap symbol
  monomerRe = /\[([A-Za-z0-9_\-,()]+)\]|(.)/g;

  /** Split sequence for single character monomers, square brackets multichar monomer names or gap symbol. */
  getSplitterAsFasta(lengthLimit) {
    const resFunc = function(seq) {
      const res = wu(seq.toString().matchAll(this.monomerRe))
        .take(lengthLimit)
        .map((ma) => {
          let mRes;
          const m = ma[0];
          if (m.length > 1)
            mRes = ma[1];
          else
            mRes = m;

          return mRes;
        }).toArray();

      return res;
    }.bind(this);
    resFunc.T = 'splitterAsFasta';
    return resFunc;
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
    const resFunc = function(seq) {
      this.helmRe.lastIndex = 0;
      const ea = this.helmRe.exec(seq.toString());
      const inSeq = ea ? ea[2] : null;

      const mmPostProcess = (mm) => {
        this.helmPp1Re.lastIndex = 0;
        const pp1M = this.helmPp1Re.exec(mm);
        if (pp1M && pp1M.length >= 2)
          return pp1M[1];
        else
          return mm;
      };

      const mmList = inSeq ? inSeq.split('.') : [];
      const mmListRes = mmList.map(mmPostProcess);
      return mmListRes;
    }.bind(this);
    resFunc.T = 'splitterAsHelm';
    return resFunc;
  }

  sample(col, n) {
    if (col.length < n)
      throw new Error('Sample source is less than n requested.');

    const idxSet = new Set();
    while (idxSet.size < n) {
      const idx = Math.floor(Math.random() * col.length);
      if (!idxSet.has(idx)) idxSet.add(idx);
    }

    return wu(idxSet).map((idx) => col.get(idx));
  }

  // -- autostart --

  //name: autostart
  //tags: autostart
  //description: Bio bootstrap
  autostart() {
    this.logger.debug('Bio: detectors.js: autostart()');

    this.autostartContextMenu();
  }

  autostartContextMenu() {
    grok.events.onContextMenu.subscribe((event) => {
      if (event.args.item && event.args.item instanceof DG.GridCell &&
        event.args.item.tableColumn && event.args.item.tableColumn.semType === DG.SEMTYPE.MACROMOLECULE
      ) {
        const contextMenu = event.args.menu;
        const cell = event.args.item.cell; // DG.Cell

        grok.functions.call('Bio:addCopyMenu', {cell: cell, menu: contextMenu})
          .catch((err) => {
            grok.shell.error(err.toString());
          });

        event.preventDefault();
        return true;
      }
    });
  }
}

async function refineSeqSplitter(col, stats, separator) {
  let invalidateRequired = false;

  const refinerList = [
    {package: 'SequenceTranslator', name: 'refineNotationProviderForHarmonizedSequence'},
  ];

  for (const refineFuncFind of refinerList) {
    try {
      const funcList = DG.Func.find(refineFuncFind);
      if (funcList.length === 0) continue;

      const funcFc = funcList[0].prepare({col: col, stats: stats, separator: separator});
      const refineRes = (await funcFc.call()).getOutputParamValue();
      invalidateRequired ||= refineRes;
    } catch (err) {
      console.error(err);
    }
  }

  if (invalidateRequired) {
    // Applying custom notation provider MUST invalidate SeqHandler
    delete col.temp[SeqTemps.seqHandler];

    for (const view of grok.shell.tableViews) {
      if (view.dataFrame === col.dataFrame)
        view.grid.invalidate();
    }
  }
}
