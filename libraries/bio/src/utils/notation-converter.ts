/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UnitsHandler} from './units-handler';
import {getSplitterForColumn, getStats, NOTATION, SeqColStats, SplitterFunc, TAGS} from './macromolecule';

/** Class for handling conversion of notation systems in Macromolecule columns */
export class NotationConverter extends UnitsHandler {
  private _splitter: SplitterFunc | null = null;

  protected get splitter(): SplitterFunc {
    if (this._splitter === null)
      this._splitter = getSplitterForColumn(this.column);
    return this._splitter;
  }

  public toFasta(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.FASTA; }

  public toSeparator(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.SEPARATOR; }

  public toHelm(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.HELM; }

  /**
   * Convert a Macromolecule column from FASTA to SEPARATOR notation
   *
   * @param {string} separator  A specific separator to be used
   * @param {string} fastaGapSymbol  Gap symbol in FASTA, '-' by default
   * @return {DG.Column}        A new column in SEPARATOR notation
   */
  private convertFastaToSeparator(separator: string, fastaGapSymbol: string | null = null): DG.Column {
    if (fastaGapSymbol === null)
      fastaGapSymbol = this.defaultGapSymbol;

    const newColumn = this.getNewColumn(NOTATION.SEPARATOR);
    // assign the values to the newly created empty column
    newColumn.init((idx: number) => {
      const fastaPolymer = this.column.get(idx);
      const fastaMonomersArray = this.splitter(fastaPolymer);
      for (let i = 0; i < fastaMonomersArray.length; i++) {
        if (fastaMonomersArray[i] === fastaGapSymbol)
          fastaMonomersArray[i] = UnitsHandler._defaultGapSymbolsDict.SEPARATOR;
      }
      return fastaMonomersArray.join(separator);
    });
    newColumn.setTag(DG.TAGS.UNITS, NOTATION.SEPARATOR);
    newColumn.setTag(TAGS.separator, separator);
    return newColumn;
  }

  /**
   * Get the wrapper strings for HELM, depending on the type of the
   * macromolecule (peptide, DNA, RNA)
   *
   * @return {string[]} Array of wrappers
   */
  private getHelmWrappers(): string[] {
    const prefix = (this.isDna()) ? 'DNA1{' :
      (this.isRna()) ? 'RNA1{' :
        (this.isPeptide()) ? 'PEPTIDE1{' :
          'Unknown'; // this case should be handled as exceptional

    if (prefix === 'Unknown')
      throw new Error('Neither peptide, nor nucleotide');

    const postfix = '}$$$';
    const leftWrapper = (this.isDna()) ? 'D(' :
      (this.isRna()) ? 'R(' : ''; // no wrapper for peptides
    const rightWrapper = (this.isDna() || this.isRna()) ? ')P' : ''; // no wrapper for peptides
    return [prefix, leftWrapper, rightWrapper, postfix];
  }

  // A helper function for converting strings to HELM
  private convertToHelmHelper(
    sourcePolymer: string,
    sourceGapSymbol: string,
    prefix: string,
    leftWrapper: string,
    rightWrapper: string,
    postfix: string
  ): string {
    const monomerArray = this.splitter(sourcePolymer);
    const monomerHelmArray: string[] = monomerArray.map((mm: string) => {
      if (mm === sourceGapSymbol)
        return UnitsHandler._defaultGapSymbolsDict.HELM;
      else
        return `${leftWrapper}${mm}${rightWrapper}`;
    });
    return `${prefix}${monomerHelmArray.join('.')}${postfix}`;
  }

  /**
   * Convert a string with SEPARATOR/FASTA notation to HELM
   *
   * @param {string} sourcePolymer  A string to be converted
   * @param {string | null} sourceGapSymbol  An optional gap symbol, set to
   * default values ('-' for FASTA and '' for SEPARATOR) unless specified
   * @return {string}  The target HELM string
   */
  public convertStringToHelm(
    sourcePolymer: string,
    sourceGapSymbol: string | null = null
  ): string {
    if (sourceGapSymbol === null)
      sourceGapSymbol = this.defaultGapSymbol;
    const [prefix, leftWrapper, rightWrapper, postfix] = this.getHelmWrappers();
    return this.convertToHelmHelper(sourcePolymer, sourceGapSymbol, prefix, leftWrapper, rightWrapper, postfix);
  }

  /**
   * Convert a column to HELM
   *
   * @param {string | null} sourceGapSymbol
   * @return {DG.Column}
   */
  private convertToHelm(sourceGapSymbol: string | null = null): DG.Column {
    if (sourceGapSymbol === null)
      sourceGapSymbol = this.defaultGapSymbol;

    const [prefix, leftWrapper, rightWrapper, postfix] = this.getHelmWrappers();

    const newColumn = this.getNewColumn(NOTATION.HELM);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const sourcePolymer = this.column.get(idx);
      return this.convertToHelmHelper(sourcePolymer, sourceGapSymbol!, prefix, leftWrapper, rightWrapper, postfix);
    });
    newColumn.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    return newColumn;
  }

  /**
   * Convert SEPARATOR column to FASTA notation
   *
   * @param {string | null} fastaGapSymbol Optional gap symbol for FASTA
   * @return {DG.Column}  Converted column
   */
  private convertSeparatorToFasta(fastaGapSymbol: string | null = null): DG.Column {
    if (fastaGapSymbol === null)
      fastaGapSymbol = UnitsHandler._defaultGapSymbolsDict.FASTA;

    const newColumn = this.getNewColumn(NOTATION.FASTA);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const separatorPolymer = this.column.get(idx);
      // items can be monomers or separators
      const separatorItemsArray = this.splitter(separatorPolymer);
      const fastaMonomersArray: string[] = [];
      for (let i = 0; i < separatorItemsArray.length; i++) {
        const item = separatorItemsArray[i];
        if (item.length === 0) {
          fastaMonomersArray.push(fastaGapSymbol!);
        } else if (item.length > 1) {
          // the case of a multi-character monomer
          const monomer = '[' + item + ']';
          fastaMonomersArray.push(monomer);
        } else {
          fastaMonomersArray.push(item);
        }
      }
      return fastaMonomersArray.join('');
    });
    newColumn.setTag(DG.TAGS.UNITS, NOTATION.FASTA);
    return newColumn;
  }

  /**
   *  Convert HELM column to FASTA/SEPARATOR
   *
   * @param {string} tgtNotation    Target notation: FASTA or SEPARATOR
   * @param {string} tgtSeparator   Optional target separator (for HELM ->
   * @param {string | null} tgtGapSymbol   Optional target gap symbol
   * SEPARATOR)
   * @return {DG.Column} Converted column
   */
  private convertHelm(tgtNotation: string, tgtSeparator?: string, tgtGapSymbol?: string): DG.Column {
    // This function must not contain calls of isDna() and isRna(), for
    // source helm columns may contain RNA, DNA and PT across different rows
    if (!tgtGapSymbol) {
      tgtGapSymbol = (this.toFasta(tgtNotation as NOTATION)) ?
        UnitsHandler._defaultGapSymbolsDict.FASTA :
        UnitsHandler._defaultGapSymbolsDict.SEPARATOR;
    }

    if (!tgtSeparator) {
      tgtSeparator = (this.toFasta(tgtNotation as NOTATION)) ? '' : this.separator;
    }

    const helmWrappersRe = /(R\(|D\(|\)|P)/g;
    const newColumn = this.getNewColumn(tgtNotation as NOTATION);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const helmPolymer = this.column.get(idx);

      // we cannot use isDna() or isRna() because source helm columns can
      // contain DNA, RNA and PT in different cells, so the corresponding
      // tags cannot be set for the whole column
      const isNucleotide = helmPolymer.startsWith('DNA') || helmPolymer.startsWith('RNA');

      // items can be monomers or helms
      const helmItemsArray = this.splitter(helmPolymer);
      const tgtMonomersArray: string[] = [];
      for (let i = 0; i < helmItemsArray.length; i++) {
        let item = helmItemsArray[i];
        if (isNucleotide)
          item = item.replace(helmWrappersRe, '');
        if (item === UnitsHandler._defaultGapSymbolsDict.HELM) {
          tgtMonomersArray.push(tgtGapSymbol!);
        } else if (this.toFasta(tgtNotation as NOTATION) && item.length > 1) {
          // the case of a multi-character monomer converted to FASTA
          const monomer = '[' + item + ']';
          tgtMonomersArray.push(monomer);
        } else {
          tgtMonomersArray.push(item);
        }
      }
      return tgtMonomersArray.join(tgtSeparator);
    });

    // TAGS.aligned is mandatory for columns of NOTATION.FASTA and NOTATION.SEPARATOR
    const splitter: SplitterFunc = getSplitterForColumn(newColumn);
    const stats: SeqColStats = getStats(newColumn, 5, splitter);
    const aligned = stats.sameLength ? 'SEQ.MSA' : 'SEQ';
    newColumn.setTag(TAGS.aligned, aligned);

    return newColumn;
  }

  private convertHelmToSeparator(): DG.Column {
    // TODO: implementatioreturn this.getNewColumn();
    return this.getNewColumn(NOTATION.SEPARATOR);
  }

  /** Dispatcher method for notation conversion
   *
   * @param {NOTATION} tgtNotation   Notation we want to convert to
   * @param {string | null} tgtSeparator   Possible separator
   * @return {DG.Column}                Converted column
   */
  public convert(tgtNotation: NOTATION, tgtSeparator: string | null = null): DG.Column {
    // possible exceptions
    if (this.notation === tgtNotation)
      throw new Error('tgt notation is invalid');
    if (this.toSeparator(tgtNotation) && tgtSeparator === null)
      throw new Error('tgt separator is not specified');

    if (this.isFasta() && this.toSeparator(tgtNotation) && tgtSeparator !== null)
      return this.convertFastaToSeparator(tgtSeparator);
    else if ((this.isFasta() || this.isSeparator()) && this.toHelm(tgtNotation))
      return this.convertToHelm();
    else if (this.isSeparator() && this.toFasta(tgtNotation))
      return this.convertSeparatorToFasta();
    else if (this.isHelm() && this.toFasta(tgtNotation)) // the case of HELM
      return this.convertHelm(tgtNotation);
    else if (this.isHelm() && this.toSeparator(tgtNotation))
      return this.convertHelm(tgtNotation, tgtSeparator!);
    else
      throw new Error('Not supported conversion ' +
        `from source notation '${this.notation}' to target notation '${tgtNotation}'.`);
  }

  public constructor(col: DG.Column) {
    super(col);
  }
}
