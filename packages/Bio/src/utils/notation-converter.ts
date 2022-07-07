import * as DG from 'datagrok-api/dg';

export const enum NOTATION {
  // these values can be changed to "user-friendly" ones later on
  FASTA = 'fasta',
  SEPARATOR = 'separator',
  HELM = 'helm'
}

export class NotationConverter {
  private _sourceColumn: DG.Column; // the column to be converted
  private _currentUnits: string; // units of the form fasta:SEQ:NT, etc.
  private _sourceNotation: NOTATION; // current notation (without :SEQ:NT, etc.)
  private _targetNotation: NOTATION;

  private get currentUnits(): string { return this._currentUnits; }
  private get sourceColumn(): DG.Column { return this._sourceColumn; }
  public get targetNotation(): NOTATION { return this._targetNotation; }

  public set targetNotation(target: NOTATION) {
    this._targetNotation = target;
  }
  public get sourceNotation(): NOTATION { return this._sourceNotation; }

  //   // these values can be changed to "user-friendly" ones later on
  //   private _fasta = 'fasta';
  //   private _separator = 'separator';
  //   private _helm = 'helm';

  public isFasta(): boolean { return this.sourceNotation === NOTATION.FASTA; }
  public isSeparator(): boolean { return this.sourceNotation === NOTATION.SEPARATOR; }
  public isHelm(): boolean { return this.sourceNotation === NOTATION.HELM; }

  private determineSourceNotation() : NOTATION {
    if (this.currentUnits.toLowerCase().startsWith('fasta'))
      return NOTATION.FASTA;
    else if (this.currentUnits.toLowerCase().startsWith('separator'))
      return NOTATION.SEPARATOR;
    else
      // TODO: handle possible exceptions
      return NOTATION.HELM;
  }

  private getNewColumn(): DG.Column {
    const col = this.sourceColumn;
    const len = col.length;
    const name = this.sourceNotation + '2' + this.targetNotation;
    const newColName = col.dataFrame.columns.getUnusedName(name);
    // dummy code
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill(name));
    newColumn.semType = 'Macromolecule';
    return newColumn;
  }

  private convertFastaToSeparator(): DG.Column {
    // TODO: implementation
    return this.getNewColumn();
  }

  private convertFastaToHelm(): DG.Column {
    // TODO: implementation
    return this.getNewColumn();
  }

  private convertSeparatorToFasta(): DG.Column {
    // TODO: implementation
    return this.getNewColumn();
  }

  private convertSeparatorToHelm(): DG.Column {
    // TODO: implementation
    return this.getNewColumn();
  }

  private convertHelmToFasta(): DG.Column {
    // TODO: implementation
    return this.getNewColumn();
  }

  private convertHelmToSeparator(): DG.Column {
    // TODO: implementatioreturn this.getNewColumn();
    return this.getNewColumn();
  }

  // TODO: write the bodies of converter methods
  public convert() : DG.Column {
    if (this.sourceNotation === this.targetNotation)
      throw new Error('Target notation is not specified');
    if (
      this.sourceNotation === NOTATION.FASTA &&
      this.targetNotation === NOTATION.SEPARATOR
    )
      return this.convertFastaToSeparator();
    else if (
      this.sourceNotation === NOTATION.FASTA &&
      this.targetNotation === NOTATION.HELM
    )
      return this.convertFastaToHelm();
    else if (
      this.sourceNotation === NOTATION.SEPARATOR &&
      this.targetNotation === NOTATION.FASTA
    )
      return this.convertSeparatorToFasta();
    else if (
      this.sourceNotation === NOTATION.SEPARATOR &&
      this.targetNotation === NOTATION.HELM
    )
      return this.convertSeparatorToHelm();
    else if (
      this.sourceNotation === NOTATION.HELM &&
      this.targetNotation === NOTATION.FASTA
    )
      return this.convertHelmToFasta();
    else
      return this.convertHelmToSeparator();
  }

  public constructor(col: DG.Column) {
    this._sourceColumn = col;
    this._currentUnits = this._sourceColumn.tags[DG.TAGS.UNITS];
    this._sourceNotation = this.determineSourceNotation();
    this._targetNotation = this.sourceNotation;
  }
}
