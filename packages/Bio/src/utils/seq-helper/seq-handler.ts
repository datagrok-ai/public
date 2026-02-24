/* eslint-disable max-lines */
/* eslint-disable @typescript-eslint/no-unused-vars */
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

/* eslint-disable max-len */
import {ALIGNMENT, ALPHABET, candidateAlphabets, getSplitterWithSeparator, NOTATION, positionSeparator, splitterAsFasta, splitterAsHelm, TAGS} from '@datagrok-libraries/bio/src/utils/macromolecule/index';
import {CandidateType, INotationProvider, ISeqConnection, ISeqSplitted, SeqColStats, SplitterFunc,} from '@datagrok-libraries/bio/src/utils/macromolecule/types';
import {detectAlphabet, detectHelmAlphabet, splitterAsFastaSimple, StringListSeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/utils';
import {mmDistanceFunctions, MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {mmDistanceFunctionType} from '@datagrok-libraries/ml/src/macromolecule-distance-functions/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_POLYMER_TYPE, HELM_WRAPPERS_REGEXP, PHOSPHATE_SYMBOL} from '@datagrok-libraries/bio/src/utils/const';
import {GAP_SYMBOL, GapOriginals} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {CellRendererBackBase, GridCellRendererTemp} from '@datagrok-libraries/bio/src/utils/cell-renderer-back-base';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {HelmType} from '@datagrok-libraries/bio/src/helm/types';
import {ConvertFunc, ISeqHandler, JoinerFunc, SeqTemps, SeqValueBase} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';

import {SeqHelper} from './seq-helper';

/** Class for handling notation units in Macromolecule columns and
 * conversion of notation systems in Macromolecule columns
 */
export class SeqHandler implements ISeqHandler {
  protected readonly _column: DG.Column; // the column to be converted
  protected readonly _units: string; // units, of the form fasta, separator
  protected readonly _notation: NOTATION; // current notation (without :SEQ:NT, etc.)
  protected readonly _defaultGapOriginal: string;
  private _notationProvider: INotationProvider | null = null;
  private _tempReadForNotProvider = false;
  protected get notationProvider(): INotationProvider | null {
    if (!this._tempReadForNotProvider) {
      this._tempReadForNotProvider = true;
      this._notationProvider = this._notationProvider ?? this._column.temp[SeqTemps.notationProvider] ?? null;
    }
    return this._notationProvider;
  };
  protected set notationProvider(v: INotationProvider | null) { this._notationProvider = v; };

  private _splitter: SplitterFunc | null = null;

  private _refinerPromise: Promise<void> = Promise.resolve();
  public get refinerPromise(): Promise<void> { return this._refinerPromise; }

  protected constructor(col: DG.Column<string>,
    private readonly seqHelper: SeqHelper,
  ) {
    if (col.type !== DG.TYPE.STRING)
      throw new Error(`Unexpected column type '${col.type}', must be '${DG.TYPE.STRING}'.`);
    this._column = col;
    const units: string | null = this._column.meta.units;
    if (!units)
      throw new Error('Units are not specified in column');
    this._units = units!;

    this._notation = this.getNotation();
    if (this.isCustom() || this.isBiln()) {
      // this.column.temp[SeqTemps.notationProvider] must be set at detector stage
      this.notationProvider = this.column.temp[SeqTemps.notationProvider] ?? null;
    }

    const defaultGapOriginal = this.isFasta() ? GapOriginals[NOTATION.FASTA] :
      this.isSeparator() ? GapOriginals[NOTATION.SEPARATOR] :
        this.isHelm() ? GapOriginals[NOTATION.HELM] :
          this.isCustom() || this.isBiln() ? (this.notationProvider?.defaultGapOriginal ?? GapOriginals[NOTATION.SEPARATOR]) :
            undefined;
    if (defaultGapOriginal == undefined)
      throw new Error(`Unexpected defaultGapOriginal for notation '${this.notation}'`);
    this._defaultGapOriginal = defaultGapOriginal;

    if (!this.column.tags.has(TAGS.aligned) || !this.column.tags.has(TAGS.alphabet) ||
      (!this.column.tags.has(TAGS.alphabetIsMultichar) && !this.isHelm() && this.alphabet === ALPHABET.UN)
    ) {
      // The following detectors and setters are to be called because the column is likely
      // as the UnitsHandler constructor was called on the column.
      if (this.isFasta())
        this.seqHelper.setUnitsToFastaColumn(this);
      else if (this.isSeparator()) {
        const separator = col.getTag(TAGS.separator);
        this.seqHelper.setUnitsToSeparatorColumn(this, separator);
      } else if (this.isHelm())
        this.seqHelper.setUnitsToHelmColumn(this);
      else if (this.isCustom() || this.isBiln())
        this.notationProvider?.setUnits(this);
      else
        throw new Error(`Unexpected units '${this.column.meta.units}'.`);
    }

    // if (!this.column.tags.has(TAGS.alphabetSize)) {
    //   if (this.isHelm())
    //     throw new Error(`For column '${this.column.name}' of notation '${this.notation}' ` +
    //       `tag '${TAGS.alphabetSize}' is mandatory.`);
    //   else if (['UN'].includes(this.alphabet))
    //     throw new Error(`For column '${this.column.name}' of alphabet '${this.alphabet}' ` +
    //       `tag '${TAGS.alphabetSize}' is mandatory.`);
    // }

    if (!this.column.tags.has(TAGS.alphabetIsMultichar)) {
      if (this.isHelm())
        this.column.setTag(TAGS.alphabetIsMultichar, 'true');
      else if (['UN'].includes(this.alphabet)) {
        throw new Error(`For column '${this.column.name}' of alphabet '${this.alphabet}' ` +
          `tag '${TAGS.alphabetIsMultichar}' is mandatory.`);
      }
    }

    this.columnVersion = this.column.version;
    // refine separator only at this stage
    if (this.isSeparator() && ((!this.isCustom() && !this.isBiln()) || !this.notationProvider) && !col.temp['seqHandlerRefined']) {
      this._refinerPromise = this.refineSeparatorNotation();
      col.temp['seqHandlerRefined'] = true;
    }
  }

  private async refineSeparatorNotation() {
    // from detectors.
    const SEQ_SAMPLE_LIMIT = 100; const SEQ_SAMPLE_LENGTH_LIMIT = 100;
    const categoriesSample = [...new Set((this.column.length < SEQ_SAMPLE_LIMIT ?
      wu.count(0).take(Math.min(SEQ_SAMPLE_LIMIT, this.column.length)).map((rowI) => this.column.get(rowI)) :
      this.column.categories.slice(0, SEQ_SAMPLE_LIMIT))
      .map((seq) => !!seq ? seq.substring(0, SEQ_SAMPLE_LENGTH_LIMIT * 5) : '')
      .filter((seq) => seq.length !== 0/* skip empty values for detector */),
    )];

    const getStats = (values: string[], minLength: number, splitter: (s: string) => string[]) => {
      const freq:{[key: string]: number} = {};
      let sameLength = true;
      let firstLength = null;

      for (const seq of values) {
        const mSeq = !!seq ? splitter(seq) : [];

        if (firstLength == null)
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
    };

    const stats = getStats(categoriesSample, 3, (s) => s.split(this.separator!));
    let invalidateRequired = false;
    const refinerList = DG.Func.find({meta: {role: 'notationRefiner'}});

    for (const refineFuncFind of refinerList) {
      try {
        const funcFc = refineFuncFind.prepare({col: this.column, stats: stats, separator: this.separator});
        const refineRes = (await funcFc.call()).getOutputParamValue();
        invalidateRequired ||= refineRes;
      } catch (err) {
        console.error(err);
      }
    }

    if (invalidateRequired) {
    // Applying custom notation provider MUST invalidate SeqHandler
      // some things might still have the old seqHandler attached, so we need to make sure they have access to notationProvider
      this._tempReadForNotProvider = false;
      delete this.column.temp[SeqTemps.seqHandler];

      this.column.fireValuesChanged();
    }
  }

  /** From detectMacromolecule */
  public static setTags(uh: SeqHandler): void {
    const units = uh.column.meta.units as NOTATION;

    if ([NOTATION.FASTA, NOTATION.SEPARATOR].includes(units)) {
      // Empty monomer alphabet is allowed, only if alphabet tag is annotated
      if (!uh.column.getTag(TAGS.alphabet) && Object.keys(uh.stats.freq).length === 0)
        throw new Error('Alphabet is empty and not annotated.');

      let aligned = uh.column.getTag(TAGS.aligned);
      if (aligned == null) {
        aligned = uh.stats.sameLength ? ALIGNMENT.SEQ_MSA : ALIGNMENT.SEQ;
        uh.column.setTag(TAGS.aligned, aligned);
      }

      let alphabet = uh.column.getTag(TAGS.alphabet);
      if (alphabet == null) {
        alphabet = detectAlphabet(uh.stats.freq, candidateAlphabets);
        uh.column.setTag(TAGS.alphabet, alphabet);
      }
      if (alphabet === ALPHABET.UN) {
        const alphabetSize = Object.keys(uh.stats.freq).length;
        const alphabetIsMultichar = Object.keys(uh.stats.freq).some((m) => m.length > 1);
        uh.column.setTag(TAGS.alphabetSize, alphabetSize.toString());
        uh.column.setTag(TAGS.alphabetIsMultichar, alphabetIsMultichar ? 'true' : 'false');
      }
    } else if (units === NOTATION.HELM) {
      let alphabet = uh.column.getTag(TAGS.alphabet);
      if (alphabet == null) {
        // const cats = uh.column.categories;
        // const splitter = uh.getSplitter();
        // const samples = Array.from(new Set(
        //   wu.count(0).take(Math.min(100, cats.length)).map((_) => Math.floor(Math.random() * cats.length)).toArray())
        // ).map((catIndex) => cats[catIndex]).filter((s) => (s?.length ?? 0) > 0).map((s) => splitter(s));
        // // splitted helm has info about polymer types
        // const polymerTypes = new Set<HELM_POLYMER_TYPE>();
        // for (const ss of samples) {
        //   ss.graphInfo?.polymerTypes
        // }
        // increase the detection threshold for candidate alphabets
        const modifiedCandidateAlphabets = candidateAlphabets.map((ca) => new CandidateType(ca.name, ca.alphabet, 0.9));

        alphabet = detectHelmAlphabet(uh.stats.freq, modifiedCandidateAlphabets, uh.defaultGapOriginal);
        uh.column.setTag(TAGS.alphabet, alphabet);
      }
    }
  }

  get column(): DG.Column { return this._column; }

  public get length(): number { return this._column.length; }

  public get units(): string { return this._units; }

  public get notation(): NOTATION { return this._notation; }

  public get defaultGapOriginal(): string { return this._defaultGapOriginal; }

  public get separator(): string | undefined {
    const separator: string | undefined = this.column.getTag(TAGS.separator) ?? undefined;
    if (this.notation === NOTATION.SEPARATOR && separator == undefined)
      throw new Error(`Separator is mandatory  for column '${this.column.name}' of notation '${this.notation}'.`);
    return separator;
  }

  public get aligned(): string {
    const aligned = this.column.getTag(TAGS.aligned);

    // TAGS.aligned is mandatory for columns of NOTATION.FASTA and NOTATION.SEPARATOR
    if (!aligned && (this.isFasta() || this.isSeparator()))
      throw new Error('Tag aligned not set');

    return aligned;
  }

  /** Alphabet name (upper case) */
  public get alphabet(): string {
    const alphabet = this.column.getTag(TAGS.alphabet);

    // TAGS.alphabet is mandatory for columns of NOTATION.FASTA and NOTATION.SEPARATOR
    if (!alphabet && (this.isFasta() || this.isSeparator()))
      throw new Error('Tag alphabet not set');

    return alphabet;
  }

  public get defaultBiotype(): HelmType {
    return this.alphabet === ALPHABET.RNA || this.alphabet === ALPHABET.DNA ? HelmTypes.NUCLEOTIDE : HelmTypes.AA;
  }

  protected get helmCompatible(): string | undefined {
    return this.column.getTag(TAGS.isHelmCompatible);
  }

  public getAlphabetSize(): number {
    if (this.notation == NOTATION.HELM || this.alphabet == ALPHABET.UN) {
      const alphabetSizeStr = this.column.getTag(TAGS.alphabetSize);
      let alphabetSize: number;
      if (alphabetSizeStr)
        alphabetSize = parseInt(alphabetSizeStr);
      else {
        // calculate alphabetSize on demand
        const stats = this.stats;
        alphabetSize = Object.keys(stats.freq).length;
      }
      return alphabetSize;
    } else {
      switch (this.alphabet) {
      case ALPHABET.PT:
        return 20;
      case ALPHABET.DNA:
      case ALPHABET.RNA:
        return 4;
      case 'NT':
        console.warn(`Unexpected alphabet 'NT'.`);
        return 4;
      default:
        throw new Error(`Unexpected alphabet '${this.alphabet}'.`);
      }
    }
  }

  public getAlphabetIsMultichar(): boolean {
    if (this.notation === NOTATION.HELM)
      return true;
    else if (this.alphabet !== ALPHABET.UN)
      return false;
    else
      return this.column.getTag(TAGS.alphabetIsMultichar) === 'true';
  }

  private cached: boolean = true;
  private _splitted: WeakRef<ISeqSplitted>[] | null = null;
  private columnVersion: number | null = null;
  // /** */
  // public get splitted(): ISeqSplitted[] {
  //   // TODO: Disable cache or invalidate on changing data
  //   if (this._splitted === null) {
  //     const splitter = this.splitter;
  //     const colLength: number = this._column.length;
  //     this._splitted = new Array(colLength);
  //     const catIdxList = this._column.getRawData();
  //     const catList: string[] = this._column.categories;
  //     for (let rowIdx: number = 0; rowIdx < colLength; rowIdx++) {
  //       const seq: string = catList[catIdxList[rowIdx]];
  //       this._splitted[rowIdx] = splitter(seq);
  //     }
  //   }
  //   return this._splitted;
  // }
  public getSplitted(rowIdx: number, limit?: number): ISeqSplitted {
    if (!this.cached || limit !== undefined) {
      const seq = this.column.get(rowIdx);
      return this.getSplitter(limit)(seq);
    } else {
      if (this.column.version !== this.columnVersion || this._splitted == null) {
        this.columnVersion = this.column.version;
        this._splitted = new Array<WeakRef<ISeqSplitted>>(this.column.length);
      }

      let resSS: ISeqSplitted | undefined = this._splitted[rowIdx] ? this._splitted[rowIdx].deref() : undefined;
      if (!resSS) {
        const seq = this.column.get(rowIdx);
        resSS = this.splitter(seq);
        this._splitted[rowIdx] = new WeakRef(resSS);
      }
      return resSS;
    }
  }

  /// Faster method to get monomers at certain position.
  /// for canonical sequences in fasta (large proteins/nucleotides)
  /// will be faster than getSplitted(rowIdx).getCanonical(posIdx)
  getMonomerAtPosition(rowIdx: number, posIdx: number, canonical: boolean): string {
    if (this.isCanonicalAlphabet && this.isFasta() && !this.getAlphabetIsMultichar()) {
      const seq = this.column.get(rowIdx) ?? '';
      const res = seq[posIdx];
      return canonical ? (res === GapOriginals[NOTATION.FASTA] ? GAP_SYMBOL : (res ?? GAP_SYMBOL)) : (res ?? '');
    }
    const mSeq: ISeqSplitted = this.getSplitted(rowIdx);
    if (posIdx < 0 || posIdx >= mSeq.length)
      return this.defaultGapOriginal;
    return canonical ? mSeq.getCanonical(posIdx) : mSeq.getOriginal(posIdx);
  }

  getMonomersAtPosition(position: number, canonical: boolean): string[] {
    const length = this.column.length;
    const res: string[] = new Array(length).fill(GAP_SYMBOL);
    if (this.isCanonicalAlphabet && this.isFasta() && !this.getAlphabetIsMultichar()) {
      const colCategories = this.column.categories;
      const colIndexes = this.column.getRawData();
      for (let i = 0; i < length; i++) {
        const seq = colCategories[colIndexes[i]] ?? '';
        if (position < seq.length) {
          const resChar = seq[position];
          res[i] = canonical ? (resChar === GapOriginals[NOTATION.FASTA] ? GAP_SYMBOL : (resChar ?? GAP_SYMBOL)) : (resChar ?? '');
        }
      }
    } else {
      for (let i = 0; i < length; i++) {
        const mSeq: ISeqSplitted = this.getSplitted(i);
        if (position >= 0 && position < mSeq.length)
          res[i] = canonical ? mSeq.getCanonical(position) : mSeq.getOriginal(position);
      }
    }
    return res;
  }

  /** Any Macromolecule can be represented on Helm format. The reverse is not always possible. */
  public getValue(rowIdx: number, options?: any): SeqValueBase {
    const seq: string = this.column.get(rowIdx);
    let resHelm: string;
    const resSeqValue = new SeqValueBase(rowIdx, this);
    return resSeqValue;
  }

  public getHelm(rowIdx: number): string {
    let resHelm: string;
    const seq = this.column.get(rowIdx);
    if (this.notation === NOTATION.HELM)
      resHelm = seq;
    else if (this.notation === NOTATION.CUSTOM || this.notation === NOTATION.BILN)
      resHelm = this.notationProvider!.getHelm(seq, {});
    else
      resHelm = this.getConverter(NOTATION.HELM)(seq);
    return resHelm;
  }

  private _stats: SeqColStats | null = null;

  public get isCanonicalAlphabet(): boolean {
    const alphabet = this.alphabet;
    return alphabet === ALPHABET.DNA || alphabet === ALPHABET.RNA || alphabet === ALPHABET.PT;
  }

  public get stats(): SeqColStats {
    if (this._stats == null) {
      const freq: { [m: string]: number } = {};
      let sameLength = true;
      let firstLength = null;

      const colLen = this.column.length;
      for (let rowIdx: number = 0; rowIdx < colLen; ++rowIdx) {
        const mSeq: ISeqSplitted = this.getSplitted(rowIdx);
        if (firstLength == null)
          firstLength = mSeq.length;
        else if (mSeq.length !== firstLength)
          sameLength = false;

        for (let posIdx = 0; posIdx < mSeq.length; ++posIdx) {
          const cm = mSeq.getCanonical(posIdx);
          if (!(cm in freq))
            freq[cm] = 0;
          freq[cm] += 1;
        }
      }
      this._stats = {freq: freq, sameLength: sameLength};
    }
    return this._stats;
  }

  private _maxLength: number | null = null;
  public get maxLength(): number {
    if (this._maxLength == null) {
      this._maxLength = this.column.length === 0 ? 0 :
        wu.count(0).take(this.column.length).map((rowIdx) => this.getSplitted(rowIdx).length).reduce((a, b) => a > b ? a : b, 0);
    }
    return this._maxLength!;
  }

  private _posList: string[] | null = null;
  public get posList(): string[] {
    if (this._posList == null) {
      const posListTxt = this.column.getTag(TAGS.positionNames);
      this._posList = posListTxt ? posListTxt.split(positionSeparator).map((p) => p.trim()) :
        wu.count(1).take(this.maxLength).map((pos) => pos.toString()).toArray();
    }
    return this._posList!;
  }

  public isFasta(): boolean { return this.notation === NOTATION.FASTA; }

  public isSeparator(): boolean { return this.notation === NOTATION.SEPARATOR || !!this.separator; }

  public isFastaOrSeparator(): boolean { return this.isFasta() || this.isSeparator(); }

  public isHelm(): boolean { return this.notation === NOTATION.HELM; }

  public isCustom(): boolean { return this.notation === NOTATION.CUSTOM; }

  public isBiln(): boolean { return this.notation === NOTATION.BILN; }

  public isRna(): boolean { return this.alphabet === ALPHABET.RNA; }

  public isDna(): boolean { return this.alphabet === ALPHABET.DNA; }

  public isPeptide(): boolean { return this.alphabet === ALPHABET.PT; }

  public isMsa(): boolean { return this.aligned ? this.aligned.toUpperCase().includes('MSA') : false; }

  public isHelmCompatible(): boolean { return this.helmCompatible === 'true'; }

  /** Checks {@link om} for being a gap
   * @param {string} om Original monomer of sequence symbol
   * @return {boolean}
   */
  public isGap(om: string): boolean {
    return !om || om === this._defaultGapOriginal;
  }

  /** Associate notation types with the corresponding units */
  /**
   * @return {NOTATION}     Notation associated with the units type
   */
  protected getNotation(): NOTATION {
    if (this.units.toLowerCase().startsWith(NOTATION.FASTA))
      return NOTATION.FASTA;
    else if (this.units.toLowerCase().startsWith(NOTATION.SEPARATOR))
      return NOTATION.SEPARATOR;
    else if (this.units.toLowerCase().startsWith(NOTATION.HELM))
      return NOTATION.HELM;
    else if (this.units.toLowerCase().startsWith(NOTATION.CUSTOM))
      return NOTATION.CUSTOM;
    else if (this.units.toLowerCase().startsWith(NOTATION.BILN))
      return NOTATION.BILN;
    else
      throw new Error(`Column '${this.column.name}' has unexpected notation '${this.units}'.`);
  }


  /**
   * Get the wrapper strings for HELM, depending on the type of the
   * macromolecule (peptide, DNA, RNA)
   *
   * @return {string[]} Array of wrappers
   */
  public getHelmWrappers(): string[] {
    const prefix = (this.isDna()) ? 'RNA1{' :
      (this.isRna() || this.isHelmCompatible()) ? 'RNA1{' : 'PEPTIDE1{';

    const postfix = '}$$$$';
    const leftWrapper = (this.isDna()) ? 'd(' :
      (this.isRna()) ? 'r(' : '';
    const rightWrapper = (this.isDna() || this.isRna()) ? ')p' : '';
    return [prefix, leftWrapper, rightWrapper, postfix];
  }

  /**
   * Create a new empty column of the specified notation type and the same
   * length as column
   *
   * @param {NOTATION} tgtNotation
   * @return {DG.Column}
   */
  protected getNewColumn(
    tgtNotation: NOTATION, tgtSeparator?: string, colName?: string, data?: string[]
  ): DG.Column<string> {
    const col = this.column;
    const name = tgtNotation.toLowerCase() + '(' + col.name + ')';
    const newColName = colName ?? col.dataFrame?.columns.getUnusedName(name) ?? name;
    const newColumn = DG.Column.fromList('string', newColName, data ?? new Array(this.column.length).fill(''));
    newColumn.semType = DG.SEMTYPE.MACROMOLECULE;
    newColumn.meta.units = tgtNotation;
    if (tgtNotation === NOTATION.SEPARATOR) {
      if (!tgtSeparator) throw new Error(`Notation \'${NOTATION.SEPARATOR}\' requires separator value.`);
      newColumn.setTag(TAGS.separator, tgtSeparator);
    }
    newColumn.setTag(DG.TAGS.CELL_RENDERER, tgtNotation === NOTATION.HELM ? 'helm' : 'sequence'); // cell.renderer

    const srcAligned = col.getTag(TAGS.aligned);
    if (srcAligned)
      newColumn.setTag(TAGS.aligned, srcAligned);

    let srcAlphabet = col.getTag(TAGS.alphabet);
    if (!srcAlphabet && this.notation === NOTATION.HELM && tgtNotation !== NOTATION.HELM)
      srcAlphabet = ALPHABET.UN;
    if (srcAlphabet != null)
      newColumn.setTag(TAGS.alphabet, srcAlphabet);

    let srcAlphabetSize: string = col.getTag(TAGS.alphabetSize);
    if (srcAlphabet != null && srcAlphabetSize)
      newColumn.setTag(TAGS.alphabetSize, srcAlphabetSize);

    const srcAlphabetIsMultichar: string = col.getTag(TAGS.alphabetIsMultichar);
    if (srcAlphabet != null && srcAlphabetIsMultichar !== undefined)
      newColumn.setTag(TAGS.alphabetIsMultichar, srcAlphabetIsMultichar);

    if (tgtNotation == NOTATION.HELM) {
      srcAlphabetSize = this.getAlphabetSize().toString();
      newColumn.setTag(TAGS.alphabetSize, srcAlphabetSize);
    }

    // if its biln, we need to set it as a separator column, later to be refined
    if (tgtNotation === NOTATION.BILN) {
      newColumn.setTag(TAGS.separator, '-');
      newColumn.meta.units = NOTATION.SEPARATOR;
    }

    return newColumn;
  }

  /** Creates a new column on data of {@link seqList} with the same tags */
  public getNewColumnFromList(name: string, seqList: string[]): DG.Column<string> {
    return this.getNewColumn(this.notation, this.separator, name, seqList);
  }

  /**
   * A helper function checking the validity of the 'units' string
   *
   * @param {string} units  the string to be validated
   * @return {boolean}
   */
  public static unitsStringIsValid(units: string): boolean {
    units = units.toLowerCase();
    const prefixes = [NOTATION.FASTA, NOTATION.SEPARATOR, NOTATION.HELM];
    const postfixes = ['rna', 'dna', 'pt'];

    const prefixCriterion = prefixes.some((p) => units.startsWith(p.toLowerCase()));
    return prefixCriterion;
  }

  /**
   * Construct a new column of semantic type MACROMOLECULE from the list of
   * specified parameters
   *
   * @param {number}    len  the length of the new column
   * @param {string}    name  the name of the new column
   * @param {string}    units  the units of the new column
   * @return {DG.Column}
   */
  public static getNewColumnFromParams(
    len: number,
    name: string,
    units: string
  ): DG.Column {
    // WARNING: in this implementation is is impossible to verify the uniqueness
    // of the new column's name
    // TODO: verify the validity of units parameter
    if (!SeqHandler.unitsStringIsValid(units))
      throw new Error('Invalid format of \'units\' parameter');
    const newColumn = DG.Column.fromList('string', name, new Array(len).fill(''));
    newColumn.semType = DG.SEMTYPE.MACROMOLECULE;
    newColumn.meta.units = units;
    return newColumn;
  }

  get splitter(): SplitterFunc {
    if (this._splitter == null)
      this._splitter = this.getSplitter();
    return this._splitter;
  }

  /** Gets function to split seq value to monomers */
  protected getSplitter(limit?: number): SplitterFunc {
    let splitter: SplitterFunc | null = null;
    splitter = this.notationProvider ? this.notationProvider.splitter : null;
    if (splitter) return splitter;

    if (this.units.toLowerCase().startsWith(NOTATION.FASTA)) {
      const alphabet: string | null = this.column.getTag(TAGS.alphabet);
      if (alphabet !== null && !this.getAlphabetIsMultichar())
        return splitterAsFastaSimple;
      else
        return splitterAsFasta;
    } else if (this.units.toLowerCase().startsWith(NOTATION.SEPARATOR))
      return getSplitterWithSeparator(this.separator!, limit);
    else if (this.units.toLowerCase().startsWith(NOTATION.HELM))
      return splitterAsHelm;
    else
      throw new Error(`Unexpected units ${this.units} .`);

    // TODO: Splitter for HELM
  }

  public split(seq: string): ISeqSplitted {
    return this.splitter(seq);
  }

  public getDistanceFunctionName(): MmDistanceFunctionsNames {
    // TODO add support for helm and separator notation
    if (!this.isFasta())
      throw new Error('Only FASTA notation is supported');
    if (this.isMsa())
      return MmDistanceFunctionsNames.HAMMING;
    switch (this.alphabet) {
    case ALPHABET.DNA:
    case ALPHABET.RNA:
      // As DNA and RNA scoring matrices are same as identity matrices(mostly),
      // we can use very fast and optimized Levenshtein distance library
      return MmDistanceFunctionsNames.LEVENSHTEIN;
    case ALPHABET.PT:
      return MmDistanceFunctionsNames.LEVENSHTEIN;
      // For default case, let's use Levenshtein distance
    default:
      return MmDistanceFunctionsNames.LEVENSHTEIN;
    }
  }

  public getDistanceFunction(): mmDistanceFunctionType {
    return mmDistanceFunctions[this.getDistanceFunctionName()]();
  }

  // checks if the separator notation is compatible with helm library
  public async checkHelmCompatibility(): Promise<boolean> {
    // check first for the column tag to avoid extra processing
    if (this.column.tags.has(TAGS.isHelmCompatible))
      return this.column.getTag(TAGS.isHelmCompatible) === 'true';

    // get the monomer lib and check against the column
    const monomerLibHelper: IMonomerLibHelper = await getMonomerLibHelper();
    const bioLib = monomerLibHelper.getMonomerLib();
    // retrieve peptides
    const peptides = bioLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE);
    // convert the peptides list to a set for faster lookup
    const peptidesSet = new Set(peptides);
    // get splitter for given separator and check if all monomers are in the lib
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    const splitterFunc = getSplitterWithSeparator(this.separator!);
    // iterate over the columns, split them and check if all monomers are in the lib
    //TODO maybe add missing threshold so that if there are not too many missing monomers
    // the column is still considered helm compatible
    const catIdxSet: Set<number> = new Set();
    const rowCount = this.column.length;
    const colRawData = this.column.getRawData();
    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
      const catI = colRawData[rowIdx];
      if (!(catI in catIdxSet)) {
        catIdxSet.add(catI);
        const seqSS = this.getSplitted(rowIdx);
        for (let posIdx = 0; posIdx < seqSS.length; ++posIdx) {
          const cm = seqSS.getCanonical(posIdx);
          if (!peptidesSet.has(cm)) {
            this.column.setTag(TAGS.isHelmCompatible, 'false');
            return false;
          }
        }
      }
    }
    this.column.setTag(TAGS.isHelmCompatible, 'true');
    return true;
  }

  // -- Notation Converter --

  public toFasta(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.FASTA; }

  public toSeparator(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.SEPARATOR; }

  public toHelm(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.HELM; }

  /**
   *  Convert HELM string to FASTA/SEPARATOR
   *
   * @param {string} srcSeq    A string to be converted
   * @param {string} tgtNotation    Target notation: FASTA or SEPARATOR
   * @param {string} tgtSeparator   Optional target separator (for HELM ->
   * @param {string | null} tgtGapOriginal   Optional target gap symbol
   * SEPARATOR)
   * @return {string} Converted string
   */
  public convertHelmToFastaSeparator(
    srcSeq: string, tgtNotation: string, tgtSeparator?: string, tgtGapOriginal?: string
  ): string {
    if (!tgtGapOriginal) {
      tgtGapOriginal = (this.toFasta(tgtNotation as NOTATION)) ?
        GapOriginals[NOTATION.FASTA] :
        GapOriginals[NOTATION.SEPARATOR];
    }

    if (!tgtSeparator)
      tgtSeparator = (this.toFasta(tgtNotation as NOTATION)) ? '' : this.separator;

    const isNucleotide = srcSeq.startsWith('RNA');
    // items can be monomers or helms
    const helmItemsArray = splitterAsHelm(srcSeq);
    const tgtMonomersArray: string[] = [];
    for (let posIdx = 0; posIdx < helmItemsArray.length; ++posIdx) {
      let om: string = helmItemsArray.getOriginal(posIdx);
      if (isNucleotide)
        om = om.replace(HELM_WRAPPERS_REGEXP, '$1');
      if (om === GapOriginals[NOTATION.HELM])
        tgtMonomersArray.push(tgtGapOriginal);
      else if (this.toFasta(tgtNotation as NOTATION) && om.length > 1) {
        // the case of a multi-character monomer converted to FASTA
        const monomer = '[' + om + ']';
        tgtMonomersArray.push(monomer);
      } else
        tgtMonomersArray.push(om);
    }
    return tgtMonomersArray.join(tgtSeparator);
  }

  /** Dispatcher method for notation conversion
   *
   * @param {NOTATION} tgtNotation   Notation we want to convert to
   * @param {string | null} tgtSeparator   Possible separator
   * @return {DG.Column}                Converted column
   */
  public convert(tgtNotation: NOTATION, tgtSeparator?: string): DG.Column<string> {
    // Get joiner from the source column units handler (this) knowing about the source sequence.
    // For example, converting DNA Helm to fasta requires removing the r(X)p decoration.
    if (!this.isCustom() || tgtNotation !== NOTATION.HELM) {
      const joiner: JoinerFunc = this.getJoiner({notation: tgtNotation, separator: tgtSeparator});
      const newColumn = this.getNewColumn(tgtNotation, tgtSeparator);
      // assign the values to the newly created empty column
      newColumn.init((rowIdx: number) => {
        const srcSS = this.getSplitted(rowIdx);
        return joiner(srcSS);
      });
      return newColumn;
    } else {
      if (!this.notationProvider || !this.notationProvider.getHelm)
        throw new Error('Notation provider with getHelm function is required for custom notation conversion.');
      const newColumn = this.getNewColumn(tgtNotation, tgtSeparator);
      newColumn.init((rowIdx: number) => {
        const srcSeq = this.column.get(rowIdx);
        return this.notationProvider!.getHelm(srcSeq, {});
      });
      return newColumn;
    }
  }

  /**
   * @param name
   * @param startIdx Start position index of the region (0-based)
   * @param endIdx   End position index of the region (0-based, inclusive)
   */
  public getRegion(startIdx: number | null, endIdx: number | null, name: string): DG.Column<string> {
    const regCol: DG.Column<string> = this.getNewColumn(this.notation, this.separator);
    regCol.name = name;

    const startIdxVal: number = startIdx ?? 0;
    const endIdxVal: number = endIdx ?? this.maxLength - 1;

    const joiner = this.getJoiner();

    const regLength = endIdxVal - startIdxVal + 1;
    const gapOM = GapOriginals[this.notation];
    regCol.init((rowI): string => {
      const seqS = this.getSplitted(rowI);
      // Custom slicing instead of array method to maintain gaps
      const regOMList: string[] = new Array<string>(regLength);
      for (let regJPos: number = 0; regJPos < regLength; ++regJPos) {
        const seqJPos = startIdxVal + regJPos;
        regOMList[regJPos] = seqJPos < seqS.length ? seqS.getOriginal(seqJPos) : gapOM;
      }
      return joiner(new StringListSeqSplitted(regOMList, gapOM));
    });

    const getRegionOfPositionNames = (str: string): string => {
      const srcPosList = str.split(',').map((p) => p.trim());
      const regPosList = new Array<string>(regLength);
      for (let regJPos: number = 0; regJPos < regLength; ++regJPos) {
        const srcJPos = startIdxVal + regJPos;
        regPosList[regJPos] = srcJPos < srcPosList.length ? srcPosList[srcJPos] : '?';
      }
      return regPosList.join(positionSeparator);
    };

    const srcPositionNamesStr = this.column.getTag(TAGS.positionNames);
    if (srcPositionNamesStr) regCol.setTag(TAGS.positionNames, getRegionOfPositionNames(srcPositionNamesStr));

    const srcPositionLabelsStr = this.column.getTag(TAGS.positionLabels);
    if (srcPositionLabelsStr) regCol.setTag(TAGS.positionLabels, getRegionOfPositionNames(srcPositionLabelsStr));

    return regCol;
  }

  private _joiner?: JoinerFunc = undefined;

  public get joiner(): JoinerFunc {
    if (!this._joiner)
      this._joiner = this.getJoiner();

    return this._joiner;
  }

  public getJoiner(opts?: { notation: NOTATION, separator?: string }): JoinerFunc {
    const notation = opts ? opts.notation : this.notation;
    const separator = opts ? opts.separator : this.separator;

    let res: JoinerFunc;
    const srcSh = this;
    switch (notation) {
    case NOTATION.FASTA: {
      res = function(srcSS: ISeqSplitted): string { return srcSh.joinToFasta(srcSS, srcSh.isHelm()); };
      break;
    }
    case NOTATION.SEPARATOR: {
      if (!separator) throw new Error(`Separator is mandatory for notation '${notation}'.`);
      res = function(srcSS: ISeqSplitted): string { return joinToSeparator(srcSS, separator, srcSh.isHelm()); };
      break;
    }
    case NOTATION.HELM: {
      const isDnaOrRna = srcSh.alphabet === ALPHABET.DNA || srcSh.alphabet === ALPHABET.RNA;
      const wrappers = srcSh.getHelmWrappers();
      res = function(srcSS: ISeqSplitted): string { return joinToHelm(srcSS, wrappers, isDnaOrRna); };
      break;
    }
    case NOTATION.BILN: {
      res = function(srcSS: ISeqSplitted): string { return joinToBiln(srcSS); };
      break;
    }
    default:
      throw new Error(`Unexpected notation '${notation}'.`);
    }

    return res;
  }

  public getConverter(tgtUnits: NOTATION, tgtSeparator: string | undefined = undefined): ConvertFunc {
    if (tgtUnits === NOTATION.SEPARATOR && !tgtSeparator)
      throw new Error(`Target separator is not specified for target units '${NOTATION.SEPARATOR}'.`);

    const srcSh = this;
    if (tgtUnits === NOTATION.FASTA)
      return function(srcSeq: string) { return srcSh.convertToFasta(srcSeq); };
    if (tgtUnits === NOTATION.HELM)
      return function(srcSeq: string) { return srcSh.convertToHelm(srcSeq); };
    else if (tgtUnits === NOTATION.SEPARATOR)
      return function(srcSeq: string) { return srcSh.convertToSeparator(srcSeq, tgtSeparator!); };
    else if (tgtUnits === NOTATION.BILN)
      return function(srcSeq: string) { return srcSh.convertToBiln(srcSeq); };
    else
      throw new Error('Unexpected target units \'' + tgtUnits + '\'.');
  }

  /** Gets a column's UnitsHandler object from temp slot or creates a new and stores it to the temp slot. */
  public static forColumn(col: DG.Column<string>, seqHelper: SeqHelper): SeqHandler {
    // TODO: Invalidate col.temp[Temps.uh] checking column's metadata
    let res = col.temp[SeqTemps.seqHandler];
    if (!res || res.columnVersion !== col.version)
      res = col.temp[SeqTemps.seqHandler] = new SeqHandler(col, seqHelper);
    return res;
  }

  // -- joiners & converters --

  private joinToFasta(seqS: ISeqSplitted, isHelm: boolean): string {
    const resMList: string[] = new Array<string>(seqS.length);
    for (let posIdx: number = 0; posIdx < seqS.length; ++posIdx) {
      const cm: string = seqS.getOriginal(posIdx);
      let om: string = seqS.getOriginal(posIdx);
      if (isHelm)
        om = om.replace(HELM_WRAPPERS_REGEXP, '$1');

      if (cm === GAP_SYMBOL)
        om = GapOriginals[NOTATION.FASTA];
      else if (cm === PHOSPHATE_SYMBOL)
        om = '';
      else if (om.length > 1)
        om = '[' + om + ']';

      resMList[posIdx] = om;
    }
    return resMList.join('');
  }

  private convertToFasta(src: string): string {
    const srcUhSplitter: SplitterFunc = this.splitter;
    const srcSS: ISeqSplitted = this.isHelm() ? this.splitterAsHelmNucl(src) : srcUhSplitter(src);
    return this.joinToFasta(srcSS, this.isHelm());
  }

  private convertToSeparator(src: string, tgtSeparator: string): string {
    const srcSS: ISeqSplitted = this.isHelm() ? this.splitterAsHelmNucl(src) : this.splitter(src);
    return joinToSeparator(srcSS, tgtSeparator, this.isHelm());
  }

  private convertToHelm(src: string): string {
    if (this.notation == NOTATION.HELM) return src;
    if (this.isCustom() && this.notationProvider?.getHelm)
      return this.notationProvider.getHelm(src, {});

    const wrappers = this.getHelmWrappers();

    const isDnaOrRna = this.isDna() || this.isRna();
    const srcSS = this.splitter(src);
    return joinToHelm(srcSS, wrappers, isDnaOrRna);
  }

  private convertToBiln(src: string): string {
    if (this.notation == NOTATION.BILN) return src;
    const srcSS = this.splitter(src);
    return joinToBiln(srcSS);
  }

  /** Splits Helm sequence adjusting nucleotides to single char symbols. (!) Removes lone phosphorus. */
  private splitterAsHelmNucl(src: string): ISeqSplitted {
    const srcMList: ISeqSplitted = this.splitter(src);
    const tgtMList: (string | null)[] = new Array<string>(srcMList.length);
    const isDna = src.startsWith('DNA');
    const isRna = src.startsWith('RNA');
    for (let posIdx: number = 0; posIdx < srcMList.length; ++posIdx) {
      let om: string | null = srcMList.getOriginal(posIdx);
      if (isDna || isRna) {
        om = om.replace(HELM_WRAPPERS_REGEXP, '$1');
        om = om === PHOSPHATE_SYMBOL ? null : om;
      }
      tgtMList[posIdx] = om ? om : null;
    }
    return new StringListSeqSplitted(tgtMList.filter((om) => !!om) as string[], GapOriginals[NOTATION.HELM]);
  }

  // Custom notation provider

  getRendererBack(gridCol: DG.GridColumn | null, tableCol: DG.Column<string>): CellRendererBackBase<string> | null {
    const temp = this.column.temp as GridCellRendererTemp<any>;
    let res = temp.rendererBack;
    if (!res)
      res = temp.rendererBack = this.notationProvider?.createCellRendererBack(gridCol, tableCol);
    return res;
  }
}

// -- joiners --

function joinToSeparator(seqS: ISeqSplitted, tgtSeparator: string, isHelm: boolean): string {
  const resMList: string[] = new Array<string>(seqS.length);
  for (let posIdx: number = 0; posIdx < seqS.length; ++posIdx) {
    const cm = seqS.getCanonical(posIdx);
    let om = seqS.getOriginal(posIdx);
    if (isHelm)
      om = om.replace(HELM_WRAPPERS_REGEXP, '$1');

    if (cm === GAP_SYMBOL)
      om = GapOriginals[NOTATION.SEPARATOR];
    else if (cm === PHOSPHATE_SYMBOL)
      om = '';
    resMList[posIdx] = om;
  }
  return resMList.join(tgtSeparator);
}

function joinToHelm(srcSS: ISeqSplitted, wrappers: string[], isDnaOrRna: boolean): string {
  if (!srcSS.graphInfo || !((srcSS.graphInfo.connections?.length ?? 0) > 0)) {
    // no graph info - linear sequence
    const [prefix, leftWrapper, rightWrapper, postfix] = wrappers;
    const resOMList: string[] = new Array<string>(srcSS.length);
    for (let posIdx: number = 0; posIdx < srcSS.length; ++posIdx) {
      const cm = srcSS.getCanonical(posIdx);
      let om: string = cm;
      if (cm === GAP_SYMBOL)
        om = GapOriginals[NOTATION.HELM];
      else {
        if (isDnaOrRna)
          om = om.replace(HELM_WRAPPERS_REGEXP, '$1');
        om = om.length === 1 ? `${leftWrapper}${om}${rightWrapper}` : `${leftWrapper}[${om}]${rightWrapper}`;
      }
      resOMList[posIdx] = om;
    }
    return `${prefix}${resOMList.join('.')}${postfix}`;
  } else {
    // there is a graph info - so we will need to be a bit tricky here
    const seqType = isDnaOrRna ? 'RNA' : 'PEPTIDE';
    const postFix = '$$$'; // three dollar signs - one is placed at the end of sequences
    const disjointSequenceIdxs = srcSS.graphInfo.disjointSeqStarts;
    const leftWrapper = wrappers[1];
    const rightWrapper = wrappers[2];
    const disjointSequences: string[] = [];
    for (let i = 0; i < disjointSequenceIdxs.length; i++) {
      const startIdx = disjointSequenceIdxs[i];
      const endIdx = i + 1 < disjointSequenceIdxs.length ? disjointSequenceIdxs[i + 1] : srcSS.length;
      const resOMList: string[] = new Array<string>(endIdx - startIdx);
      for (let posIdx = startIdx; posIdx < endIdx; ++posIdx) {
        const cm = srcSS.getCanonical(posIdx);
        let om: string = cm;
        if (cm === GAP_SYMBOL)
          om = GapOriginals[NOTATION.HELM];
        else {
          if (isDnaOrRna)
            om = om.replace(HELM_WRAPPERS_REGEXP, '$1');
          om = om.length === 1 ? `${leftWrapper}${om}${rightWrapper}` : `${leftWrapper}[${om}]${rightWrapper}`;
        }
        resOMList[posIdx - startIdx] = om;
      }
      disjointSequences.push(`${seqType}${i + 1}{${resOMList.join('.')}}`);
    }
    // PEPTIDE2,PEPTIDE2,16:R2-1:R1|PEPTIDE3,PEPTIDE3,16:R2-1:R1|PEPTIDE3,PEPTIDE2,10:R3-1:R3|PEPTIDE1,PEPTIDE2,1:R2-9:R3$$$V2.0
    const sequencePart = disjointSequences.join('|');

    const sequenceConnections = srcSS.graphInfo.connections.map((conn) => {
      return `${seqType}${conn.seqIndex1 + 1},${seqType}${conn.seqIndex2 + 1},${conn.monomerIndex1 + 1}:R${conn.rGroup1}-${conn.monomerIndex2 + 1}:R${conn.rGroup2}`;
    }).join('|');

    return `${sequencePart}$${sequenceConnections}${postFix}V2.0`;
  }
}

function joinToBiln(srcSS: ISeqSplitted): string {
  const needsSquareBrackets = (cm: string | null) => {
    return cm && (cm.includes('-') || cm.includes('*') || cm.includes('[R'));
  };

  if (!srcSS.graphInfo || !((srcSS.graphInfo.connections?.length ?? 0) > 0)) {
    const resOMList: string[] = new Array<string>(srcSS.length);
    for (let posIdx: number = 0; posIdx < srcSS.length; ++posIdx) {
      const canonical = srcSS.getCanonical(posIdx);
      if (needsSquareBrackets(canonical)) // Biln uses '-' as a separator, need to enclose in []. also there might be smiles in there, where Rs are represented as '*' or R
        resOMList[posIdx] = `[${canonical}]`;
      else
        resOMList[posIdx] = canonical;
    }
    return resOMList.join('-'); // Biln uses '-' as a separator
  } else { // conversion happens only if there is a graph info
    const disjointSequenceIdxs = srcSS.graphInfo.disjointSeqStarts;
    const allSeqParts = new Array<string>(srcSS.length);
    for (let posIdx = 0; posIdx < srcSS.length; ++posIdx) {
      const canonical = srcSS.getCanonical(posIdx);
      // allSeqParts[posIdx] = srcSS.getCanonical(posIdx);
      if (needsSquareBrackets(canonical)) // Biln uses '-' as a separator, need to enclose in []
        allSeqParts[posIdx] = `[${canonical}]`;
      else
        allSeqParts[posIdx] = canonical;
    }
    for (let i = 0; i < srcSS.graphInfo.connections.length; i++) {
      const conn: ISeqConnection = srcSS.graphInfo.connections[i];
      const conId = `${i + 1}`;
      const seq1Idx = conn.seqIndex1;
      const seq2Idx = conn.seqIndex2;
      const monomer1Idx = disjointSequenceIdxs[seq1Idx] + conn.monomerIndex1;
      const monomer2Idx = disjointSequenceIdxs[seq2Idx] + conn.monomerIndex2;
      const seqPart1 = `${allSeqParts[monomer1Idx]}(${conId},${conn.rGroup1})`;
      const seqPart2 = `${allSeqParts[monomer2Idx]}(${conId},${conn.rGroup2})`;
      allSeqParts[monomer1Idx] = seqPart1;
      allSeqParts[monomer2Idx] = seqPart2;
    }

    const disjointParts = disjointSequenceIdxs.map((startIdx, i) => {
      const endIdx = i + 1 < disjointSequenceIdxs.length ? disjointSequenceIdxs[i + 1] : srcSS.length;
      return allSeqParts.slice(startIdx, endIdx).join('-');
    });
    return disjointParts.join('.'); // Biln uses '-' as a separator and '.' between disjoint sequences
  }
}

