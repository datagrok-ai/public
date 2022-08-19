import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {WebLogo, SeqColStats} from '../viewers/web-logo';

/** enum type to simplify setting "user-friendly" notation if necessary */
export const enum NOTATION {
  FASTA = 'FASTA',
  SEPARATOR = 'SEPARATOR',
  HELM = 'HELM'
}

/** Class for handling notation units in Macromolecule columns */
export class UnitsHandler {
  protected readonly _column: DG.Column; // the column to be converted
  protected _units: string; // units, of the form fasta, separator
  protected _notation: NOTATION; // current notation (without :SEQ:NT, etc.)
  protected _defaultGapSymbol: string;
  protected _defaultGapSymbolsDict = {
    HELM: '*',
    SEPARATOR: '',
    FASTA: '-',
  };

  public static readonly PeptideFastaAlphabet = new Set([
    'G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
    'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T',
  ]);
  public static readonly DnaFastaAlphabet = new Set(['A', 'C', 'G', 'T']);
  public static readonly RnaFastaAlphabet = new Set(['A', 'C', 'G', 'U']);

  public static setUnitsToFastaColumn(col: DG.Column) {
    if (col.semType !== DG.SEMTYPE.MACROMOLECULE)
      throw new Error('Fasta column must be MACROMOLECULE');

    const stats: SeqColStats = WebLogo.getStats(col, 5, WebLogo.splitterAsFasta);
    const seqType = stats.sameLength ? 'SEQ.MSA' : 'SEQ';

    const alphabetCandidates: [string, Set<string>][] = [
      ['PT', UnitsHandler.PeptideFastaAlphabet],
      ['DNA', UnitsHandler.DnaFastaAlphabet],
      ['RNA', UnitsHandler.RnaFastaAlphabet],
    ];

    // Calculate likelihoods for alphabet_candidates
    const alphabetCandidatesSim: number[] = alphabetCandidates.map(
      (c) => WebLogo.getAlphabetSimilarity(stats.freq, c[1]));
    const maxCos = Math.max(...alphabetCandidatesSim);
    const alphabet = maxCos > 0.65 ? alphabetCandidates[alphabetCandidatesSim.indexOf(maxCos)][0] : 'UN';
    
    const units: string = 'fasta';
    col.setTag(DG.TAGS.UNITS, units);
    col.setTag('aligned', seqType);
    col.setTag('alphabet', alphabet);
  }

  protected get units(): string { return this._units; }

  protected get column(): DG.Column { return this._column; }

  public get notation(): NOTATION { return this._notation; }

  public get defaultGapSymbol(): string { return this._defaultGapSymbol; }

  public get separator(): string {
    const separator = this.column.getTag('separator');
    if (separator !== null)
      return separator;
    else
      throw new Error('Separator not set');
  }

  public get aligned(): string {
    const aligned = this.column.getTag('aligned');
    if (aligned !== null) {
      return aligned;
    } else {
      throw new Error('Tag aligned not set');
    }
  }

  public get alphabet(): string {
    const alphabet = this.column.getTag('alphabet');
    if (alphabet !== null) {
      return alphabet;
    } else {
      throw new Error('Tag alphabet not set');
    }
  }

  public isFasta(): boolean { return this.notation === NOTATION.FASTA; }

  public isSeparator(): boolean { return this.notation === NOTATION.SEPARATOR; }

  public isHelm(): boolean { return this.notation === NOTATION.HELM; }

  public isRna(): boolean { return this.alphabet.toLowerCase() === 'rna'; }

  public isDna(): boolean { return this.alphabet.toLowerCase() === 'dna'; }

  public isPeptide(): boolean { return this.alphabet.toLowerCase() === 'pt'; }

  /** Associate notation types with the corresponding units */
  /**
   * @return {NOTATION}     Notation associated with the units type
   */
  protected getNotation(): NOTATION {
    if (this.units.toLowerCase().startsWith('fasta'))
      return NOTATION.FASTA;
    else if (this.units.toLowerCase().startsWith('separator'))
      return NOTATION.SEPARATOR;
    else if (this.units.toLowerCase().startsWith('helm'))
      return NOTATION.HELM;
    else
      throw new Error('The column has units that do not correspond to any notation');
  }

  /**
   * Create a new empty column of the specified notation type and the same
   * length as column
   *
   * @param {NOTATION} targetNotation
   * @return {DG.Column}
   */
  protected getNewColumn(targetNotation: NOTATION): DG.Column {
    const col = this.column;
    const len = col.length;
    const name = targetNotation.toLowerCase() + '(' + col.name + ')';
    const newColName = col.dataFrame.columns.getUnusedName(name);
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill(''));
    newColumn.semType = DG.SEMTYPE.MACROMOLECULE;
    newColumn.setTag(
      DG.TAGS.UNITS,
      this.units.replace(
        this.notation.toLowerCase().toString(),
        targetNotation.toLowerCase().toString()
      )
    );
    newColumn.setTag(DG.TAGS.CELL_RENDERER, 'Macromolecule');

    return newColumn;
  }

  /**
   * Create a new empty column using templateCol as a template
   *
   * @param {DG.Column} templateCol  the properties and units of this column are used as a
   * template to build the new one
   * @return {DG.Column}
   */
  public static getNewColumn(templateCol: DG.Column): DG.Column {
    const col: UnitsHandler = new UnitsHandler(templateCol);
    const targetNotation = col.notation;
    return col.getNewColumn(targetNotation);
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
    if (!UnitsHandler.unitsStringIsValid(units))
      throw new Error('Invalid format of \'units\' parameter');
    const newColumn = DG.Column.fromList('string', name, new Array(len).fill(''));
    newColumn.semType = DG.SEMTYPE.MACROMOLECULE;
    newColumn.setTag(DG.TAGS.UNITS, units);
    return newColumn;
  }

  public constructor(col: DG.Column) {
    this._column = col;
    const units = this._column.tags[DG.TAGS.UNITS];
    if (units !== null)
      this._units = units;
    else
      throw new Error('Units are not specified in column');
    this._notation = this.getNotation();
    this._defaultGapSymbol = (this.isFasta()) ? this._defaultGapSymbolsDict.FASTA :
      (this.isHelm()) ? this._defaultGapSymbolsDict.HELM :
        this._defaultGapSymbolsDict.SEPARATOR;
  }
}
