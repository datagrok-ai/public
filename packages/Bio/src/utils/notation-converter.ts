import * as DG from 'datagrok-api/dg';

// export const enum NOTATION {
//   // these values can be changed to "user-friendly" ones later on
//   FASTA = 'fasta',
//   SEPARATOR = 'separator',
//   HELM = 'helm'
// }

export class NotationConverter {
  private _sourceColumn: DG.Column; // the column to be converted
  private _currentUnits: string; // units of the form fasta:SEQ:NT, etc.
  private _sourceNotation: string; // current notation (without :SEQ:NT, etc.)
  private _targetNotation: string;

  private get sourceColumn(): DG.Column { return this._sourceColumn; }
  private get currentUnits(): string { return this._currentUnits; }
  private get sourceNotation(): string { return this._sourceNotation; }
  private get targetNotation(): string { return this._targetNotation; }

  // these values can be changed to "user-friendly" ones later on
  private _fasta = 'fasta';
  private _separator = 'separator';
  private _helm = 'helm';

  public isFasta(): boolean { return this.sourceNotation == this._fasta; }
  public isSeparator(): boolean { return this.sourceNotation == this._separator; }
  public isHelm(): boolean { return this.sourceNotation == this._helm; }

  private determineSourceNotation() : string {
    if (this.currentUnits.toLowerCase().startsWith('fasta'))
      return 'fasta';
    else if (this.currentUnits.toLowerCase().startsWith('separator'))
      return 'separator';
    else
      // TODO: handle possible exceptions
      return 'HELM';
  }

  private convertFastaToSeparator(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('fasta2sep'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertFastaToHelm(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('fasta2helm'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertSeparatorToFasta(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('sep2fasta'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertSeparatorToHelm(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('sep2helm'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertHelmToFasta(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('helm2fasta'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertHelmToSeparator(): DG.Column {
    // TODO: implementation
    const len = this.sourceColumn.length;
    const newColName = 'converted';
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill('helm2sep'));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  // TODO: write the bodies of converter methods
  public convert() : DG.Column {
    if (
      this.sourceNotation == this._fasta &&
      this.targetNotation == this._separator
    )
      return this.convertFastaToSeparator();
    else if (
      this.sourceNotation == this._fasta &&
      this.targetNotation == this._helm
    )
      return this.convertFastaToHelm();
    else if (
      this.sourceNotation == this._separator &&
      this.targetNotation == this._fasta
    )
      return this.convertSeparatorToFasta();
    else if (
      this.sourceNotation == this._separator &&
      this.targetNotation == this._helm
    )
      return this.convertSeparatorToHelm();
    else if (
      this.sourceNotation == this._helm &&
      this.targetNotation == this._fasta
    )
      return this.convertHelmToFasta();
    else
      return this.convertHelmToSeparator();
  }

  public constructor(col: DG.Column, target: string) {
    this._sourceColumn = col;
    this._currentUnits = this._sourceColumn.tags[DG.TAGS.UNITS];
    this._sourceNotation = this.determineSourceNotation();
    this._targetNotation = target;
  }
}
