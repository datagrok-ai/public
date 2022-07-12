import * as DG from 'datagrok-api/dg';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';

/** enum type to simplify setting "user-friendly" notation if necessary */
export const enum NOTATION {
  FASTA = 'FASTA',
  SEPARATOR = 'SEPARATOR',
  HELM = 'HELM'
}

/** Class for handling conversion of notation systems in Macromolecule columns */
export class NotationConverter {
  private _sourceColumn: DG.Column; // the column to be converted
  private _sourceUnits: string; // units, of the form fasta:SEQ:NT, etc.
  private _sourceNotation: NOTATION; // current notation (without :SEQ:NT, etc.)
  private _defaultGapSymbol: string;
  private _defaultGapSymbolsDict = {
    helm: '*',
    separator: '',
    fasta: '-',
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
    const name = targetNotation + '(' + col.name + ')';
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
        // TODO: replace by the enumeration value
        'Macromolecule');
    }
    return newColumn;
  }

  /**
   * Convert a Macromolecule column from FASTA to SEPARATOR notation
   *
   * @param {string} separator  A specific separator to be used
   * @param {string} gapSymbol  Gap symbol in FASTA, '-' by default
   * @return {DG.Column}        A new column in SEPARATOR notation
   */
  private convertFastaToSeparator(separator: string, gapSymbol: string = '-'): DG.Column {
    // a function splitting FASTA sequence into an array of monomers:
    const splitterAsFasta = WebLogo.splitterAsFasta;

    const newColumn = this.getNewColumn(NOTATION.SEPARATOR);
    // assign the values to the newly created empty column
    newColumn.init((idx: number) => {
      const fastaPolymer = this.sourceColumn.get(idx);
      const fastaMonomersArray = splitterAsFasta(fastaPolymer);
      for (let i = 0; i < fastaMonomersArray.length; i++) {
        if (fastaMonomersArray[i] === gapSymbol)
          fastaMonomersArray[i] = '';
      }
      return fastaMonomersArray.join(separator);
    });
    return newColumn;
  }

  private convertToHelm(sourceGapSymbol: string | null = null) {
    if (sourceGapSymbol === null)
      sourceGapSymbol = this.defaultGapSymbol;
    // A function splitting a sequence into an array of monomers according to
    // its notation
    const splitter = WebLogo.getSplitterForColumn(this.sourceColumn);

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

    const newColumn = this.getNewColumn(NOTATION.HELM);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const sourcePolymer = this.sourceColumn.get(idx);
      const sourceMonomersArray = splitter(sourcePolymer);
      const helmArray = [prefix];
      let firstIteration = true;
      for (let i = 0; i < sourceMonomersArray.length; i++) {
        const dot = firstIteration ? '' : '.';
        let token = sourceMonomersArray[i];
        if (token === sourceGapSymbol)
          token = this._defaultGapSymbolsDict.helm;
        const item = [dot, leftWrapper, token, rightWrapper];
        helmArray.push(item.join(''));
        firstIteration = false;
      }
      helmArray.push(postfix);
      return helmArray.join('');
    });
    return newColumn;
  }

  private handleSeparatorItemForFasta(
    idx: number,
    separatorItemsArray: string[],
    separator: string,
    gapSymbol: string,
    fastaMonomersArray: string[]
  ): void {
    const item = separatorItemsArray[idx];
    if (item.length > 1) {
      // the case of a multi-character monomer
      const monomer = '[' + item + ']';
      fastaMonomersArray.push(monomer);
    }
    if (item === separator) {
      if (idx !== 0 && separatorItemsArray[idx - 1] === separator)
        fastaMonomersArray.push(gapSymbol);
    }
  }

  private convertSeparatorToFasta(
    separator: string | null = null,
    gapSymbol: string = '-'
  ): DG.Column {
    // TODO: implementation
    // * similarly to fasta2separator, divide string into monomers
    // * adjacent separators is a gap (symbol to be specified)
    // * the monomers MUST be single-character onles, otherwise forbid
    // * NO, they can be multi-characters
    // conversion
    // * consider automatic determining the separator

    if (separator === null)
      separator = this.separator;

    // a function splitting FASTA sequence into an array of monomers
    //const splitterAsSeparator = WebLogo.getSplitterWithSeparator(separator);
    const splitter = WebLogo.getSplitterForColumn(this.sourceColumn);

    const newColumn = this.getNewColumn(NOTATION.FASTA);
    // assign the values to the empty column
    newColumn.init((idx: number) => {
      const separatorPolymer = this.sourceColumn.get(idx);
      // items can be monomers or separators
      const separatorItemsArray = splitter(separatorPolymer);
      const fastaMonomersArray: string[] = [];
      for (let i = 0; i < separatorItemsArray.length; i++) {
        const item = separatorItemsArray[i];
        if (item.length === 0) {
          fastaMonomersArray.push(gapSymbol);
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

  private convertHelmToFasta(): DG.Column {
    // TODO: implementation
    return this.getNewColumn(NOTATION.FASTA);
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
  public convert(targetNotation: NOTATION, tgtSeparator: string | null = null): DG.Column {
    // possible exceptions
    if (this.sourceNotation === targetNotation)
      throw new Error('Target notation is invalid');
    if (this.toSeparator(targetNotation) && tgtSeparator === null)
      throw new Error('Target separator is not specified');

    if (this.isFasta() && this.toSeparator(targetNotation) && tgtSeparator !== null)
      return this.convertFastaToSeparator(tgtSeparator);
    else if ((this.isFasta() || this.isSeparator()) && this.toHelm(targetNotation))
      return this.convertToHelm();
    else if (this.isSeparator() && this.toFasta(targetNotation))
      return this.convertSeparatorToFasta(tgtSeparator!);
    else if (this.isHelm() && this.toFasta(targetNotation))
      return this.convertHelmToFasta();
    else
      return this.convertHelmToSeparator();
  }

  public constructor(col: DG.Column) {
    this._sourceColumn = col;
    const units = this._sourceColumn.tags[DG.TAGS.UNITS];
    if (units !== null)
      this._sourceUnits = units;
    else
      throw new Error('Units are not specified in column');
    this._sourceNotation = this.getSourceNotation();
    this._defaultGapSymbol = (this.isFasta()) ? this._defaultGapSymbolsDict.fasta :
      (this.isHelm()) ? this._defaultGapSymbolsDict.helm :
        this._defaultGapSymbolsDict.separator;
  }
}
