import * as DG from 'datagrok-api/dg';
import {SplitterFunc, WebLogo} from '../viewers/web-logo';

/** enum type to simplify setting "user-friendly" notation if necessary */
export const enum NOTATION {
  FASTA = 'FASTA',
  SEPARATOR = 'SEPARATOR',
  HELM = 'HELM'
}

/** Class for handling conversion of notation systems in Macromolecule columns */
export class NotationConverter {
  private readonly _sourceColumn: DG.Column; // the column to be converted
  private _sourceUnits: string; // units, of the form fasta:SEQ:NT, etc.
  private _sourceNotation: NOTATION; // current notation (without :SEQ:NT, etc.)
  private _defaultGapSymbol: string;
  private _defaultGapSymbolsDict = {
    HELM: '*',
    SEPARATOR: '',
    FASTA: '-',
  };

  private _splitter: SplitterFunc | null = null;
  protected get splitter(): SplitterFunc {
    if (this._splitter === null)
      this._splitter = WebLogo.getSplitterForColumn(this._sourceColumn);
    return this._splitter;
  };


  private get sourceUnits(): string { return this._sourceUnits; }

  private get sourceColumn(): DG.Column { return this._sourceColumn; }

  public get sourceNotation(): NOTATION { return this._sourceNotation; }

  public get defaultGapSymbol(): string { return this._defaultGapSymbol; }

  public get separator(): string {
    const separator = this.sourceColumn.getTag('separator');
    if (separator !== null)
      return separator;
    else
      throw new Error('Separator not set');
  }

  public isFasta(): boolean { return this.sourceNotation === NOTATION.FASTA; }

  public isSeparator(): boolean { return this.sourceNotation === NOTATION.SEPARATOR; }

  public isHelm(): boolean { return this.sourceNotation === NOTATION.HELM; }

  public toFasta(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.FASTA; }

  public toSeparator(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.SEPARATOR; }

  public toHelm(targetNotation: NOTATION): boolean { return targetNotation === NOTATION.HELM; }

  public isRna(): boolean { return this.sourceUnits.toLowerCase().endsWith('rna'); }

  public isDna(): boolean { return this.sourceUnits.toLowerCase().endsWith('dna'); }

  public isPeptide(): boolean { return this.sourceUnits.toLowerCase().endsWith('pt'); }

  /** Associate notation types with the corresponding units */
  /**
   * @return {NOTATION}     Notation associated with the units type
   */
  private getSourceNotation(): NOTATION {
    if (this.sourceUnits.toLowerCase().startsWith('fasta'))
      return NOTATION.FASTA;
    else if (this.sourceUnits.toLowerCase().startsWith('separator'))
      return NOTATION.SEPARATOR;
    else if (this.sourceUnits.toLowerCase().startsWith('helm'))
      return NOTATION.HELM;
    else
      throw new Error('The column has units that do not correspond to any notation');
  }

  /**
   * Create a new empty column of the specified notation type and the same
   * length as sourceColumn
   *
   * @param {NOTATION} targetNotation
   * @return {DG.Column}
   */
  private getNewColumn(targetNotation: NOTATION): DG.Column {
    const col = this.sourceColumn;
    const len = col.length;
    const name = targetNotation.toLowerCase() + '(' + col.name + ')';
    const newColName = col.dataFrame.columns.getUnusedName(name);
    // dummy code
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill(''));
    newColumn.semType = DG.SEMTYPE.MACROMOLECULE;
    newColumn.setTag(
      DG.TAGS.UNITS,
      this.sourceUnits.replace(
        this.sourceNotation.toLowerCase().toString(),
        targetNotation.toLowerCase().toString()
      )
    );
    // TODO: specify cell renderers for all cases
    if (this.toFasta(targetNotation)) {
      newColumn.setTag(
        DG.TAGS.CELL_RENDERER,
        'Macromolecule');
    }
    return newColumn;
  }

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
      const fastaPolymer = this.sourceColumn.get(idx);
      const fastaMonomersArray = this.splitter(fastaPolymer);
      for (let i = 0; i < fastaMonomersArray.length; i++) {
        if (fastaMonomersArray[i] === fastaGapSymbol)
          fastaMonomersArray[i] = this._defaultGapSymbolsDict.SEPARATOR;
      }
      return fastaMonomersArray.join(separator);
    });
    newColumn.setTag('separator', separator);
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
        return this._defaultGapSymbolsDict.HELM;
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
  ) : string {
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
      const sourcePolymer = this.sourceColumn.get(idx);
      return this.convertToHelmHelper(sourcePolymer, sourceGapSymbol!, prefix, leftWrapper, rightWrapper, postfix);
    });
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
      fastaGapSymbol = this._defaultGapSymbolsDict.FASTA;

    const newColumn = this.getNewColumn(NOTATION.FASTA);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const separatorPolymer = this.sourceColumn.get(idx);
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
  private convertHelm(
    tgtNotation: string,
    tgtSeparator: string = '',
    tgtGapSymbol: string | null = null
  ): DG.Column {
    if (tgtGapSymbol === null) {
      tgtGapSymbol = (this.toFasta(tgtNotation as NOTATION)) ?
        this._defaultGapSymbolsDict.FASTA :
        this._defaultGapSymbolsDict.SEPARATOR;
    }

    if (this.toSeparator(tgtNotation as NOTATION) && tgtSeparator === '')
      tgtSeparator = this.separator;

    const newColumn = this.getNewColumn(tgtNotation as NOTATION);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const helmPolymer = this.sourceColumn.get(idx);
      // items can be monomers or helms
      const helmItemsArray = this.splitter(helmPolymer);
      const tgtMonomersArray: string[] = [];
      for (let i = 0; i < helmItemsArray.length; i++) {
        const item = helmItemsArray[i];
        if (item === this._defaultGapSymbolsDict.HELM) {
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
    return newColumn;
  }

  private convertHelmToSeparator(): DG.Column {
    // TODO: implementatioreturn this.getNewColumn();
    return this.getNewColumn(NOTATION.SEPARATOR);
  }

  /** Dispatcher method for notation conversion
   *
   * @param {NOTATION} targetNotation   Notation we want to convert to
   * @param {string | null} tgtSeparator   Possible separator
   * @return {DG.Column}                Converted column
   */
  public convert(tgtNotation: NOTATION, tgtSeparator: string | null = null): DG.Column {
    // possible exceptions
    if (this.sourceNotation === tgtNotation)
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
    else // this.isHelm() && this.toSeparator(tgtNotation)
      return this.convertHelm(tgtNotation, tgtSeparator!);
  }

  public constructor(col: DG.Column) {
    this._sourceColumn = col;
    const units = this._sourceColumn.tags[DG.TAGS.UNITS];
    if (units !== null)
      this._sourceUnits = units;
    else
      throw new Error('Units are not specified in column');
    this._sourceNotation = this.getSourceNotation();
    this._defaultGapSymbol = (this.isFasta()) ? this._defaultGapSymbolsDict.FASTA :
      (this.isHelm()) ? this._defaultGapSymbolsDict.HELM :
        this._defaultGapSymbolsDict.SEPARATOR;
  }
}
