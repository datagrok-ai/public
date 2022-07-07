import * as DG from 'datagrok-api/dg';

export class NotationConverter {
  private _sourceColumn: DG.Column; // the column we need to translate
  private _currentUnits: string; // current units of a column
  private _sourceNotation: string; // current notation (without .SEQ:NT, etc.)
  private _targetNotation: string;

  // these values can be changed to "user-friendly" ones later on
  private _fasta = 'fastaNotation';
  private _separator = 'separatorNotation';
  private _helm = 'HELMNotation';

  public currentUnits(): string { return this._currentUnits; }
  private sourceNotation(): string { return this._sourceNotation; }

  // notation getters
  public isFasta(): boolean { return this._currentUnits == this._fasta; }
  public isSeparator(): boolean { return this._currentUnits == this._separator; }
  public isHelm(): boolean { return this._currentUnits == this._helm; }

  // determine the source notation used in the column
  private determineSourceNotation() : string {
    if (this._currentUnits.startsWith('fasta'))
      return 'fastaNotation';
    else if (this._currentUnits.startsWith('separator'))
      return 'separatorNotatoin';
    else
      // TODO: handle possible exceptions
      return 'HELMNotation';
  }

  public constructor(col: DG.Column) {
    this._sourceColumn = col;
    this._currentUnits = this._sourceColumn.tags[DG.TAGS.UNITS];
    this._sourceNotation = this.determineSourceNotation();
    this._targetNotation = ''; // to be set separately
  }
}
