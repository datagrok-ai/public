import * as DG from 'datagrok-api/dg';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';

/** enum type to simplify setting "user-friendly" notation if necessary */
export const enum NOTATION {
  FASTA = 'fasta',
  SEPARATOR = 'separator',
  HELM = 'helm'
}

/** Class for handling conversion of notation systems in Macromolecule columns */
export class NotationConverter {
  private _sourceColumn: DG.Column; // the column to be converted
  private _sourceUnits: string; // units, of the form fasta:SEQ:NT, etc.
  private _sourceNotation: NOTATION; // current notation (without :SEQ:NT, etc.)
  private _targetNotation: NOTATION;

  private get sourceUnits(): string { return this._sourceUnits; }
  private get sourceColumn(): DG.Column { return this._sourceColumn; }

  public get targetNotation(): NOTATION { return this._targetNotation; }
  public set targetNotation(target: NOTATION) { this._targetNotation = target; }
  public get sourceNotation(): NOTATION { return this._sourceNotation; }

  public isFasta(): boolean { return this.sourceNotation === NOTATION.FASTA; }
  public isSeparator(): boolean { return this.sourceNotation === NOTATION.SEPARATOR; }
  public isHelm(): boolean { return this.sourceNotation === NOTATION.HELM; }

  public toFasta(): boolean { return this.targetNotation === NOTATION.FASTA; }
  public toSeparator(): boolean { return this.targetNotation === NOTATION.SEPARATOR; }
  public toHelm(): boolean { return this.targetNotation === NOTATION.HELM; }

  // TODO: isRna
  public isRna(): boolean { return this.sourceUnits.toLowerCase().endsWith('nt'); }
  public isPeptide(): boolean { return this.sourceUnits.toLowerCase().endsWith('pt'); }

  /** Associate notation types with the corresponding units */
  /**
   * @return {NOTATION} notation associated with the units type
   */
  private determineSourceNotation() : NOTATION {
    if (this.sourceUnits.toLowerCase().startsWith('fasta'))
      return NOTATION.FASTA;
    else if (this.sourceUnits.toLowerCase().startsWith('separator'))
      return NOTATION.SEPARATOR;
    else
      // TODO: handle possible exceptions
      return NOTATION.HELM;
  }

  // TODO: write doc
  private getNewColumn(): DG.Column {
    const col = this.sourceColumn;
    const len = col.length;
    const name = this.sourceNotation + '2' + this.targetNotation;
    const newColName = col.dataFrame.columns.getUnusedName(name);
    // dummy code
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill(''));
    newColumn.semType = 'Macromolecule';
    // TODO: determine all the qualifiers (units, ...), perhaps, using detectors
    return newColumn;
  }

  // TODO: write doc
  private convertFastaToSeparator(): DG.Column {
    // TODO: implementation
    // * specify separator
    // * fasta gap symbol should NOT be considered '-' only, but set as a parameter
    // * in fasta every position is a monomer (no multi-char monomers), call splitToMonomers() method
    // * in the resulting jagged array, every gap symbol is to be replaced by
    // the empty string, while the monomers, to be separated by the separator
    // (specified as a parameter)
    // On splitToMonomers(): /libraries/bio/src/viewers/WebLogo --> getSplitter

    const separator = '.'; // to be specified as an argument
    const gapSymbol = '-'; // to be specified as an argument
    const splitterAsFasta = WebLogo.splitterAsFasta;
    const newColumn = this.getNewColumn();
    newColumn.init((idx: number) => {
      const sourcePolymer = this.sourceColumn.get(idx);
      const monomersArray = splitterAsFasta(sourcePolymer);
      for (let i = 0; i < monomersArray.length; i++) {
        if (monomersArray[i] === gapSymbol)
          monomersArray[i] = '';
      }
      return monomersArray.join(separator);
    });
    return newColumn;
  }

  private wrapRnaNucleotideToHelm(monomer: string) {
    
  }

  private convertFastaToHelm(): DG.Column {
    const gapSymbol = '-'; // to be specified as an argument
    const splitterAsFasta = WebLogo.splitterAsFasta;
    const newColumn = this.getNewColumn();
    newColumn.init((idx: number) => {
      const sourcePolymer = this.sourceColumn.get(idx);
      const monomersArray = splitterAsFasta(sourcePolymer);
      for (let i = 0; i < monomersArray.length; i++) {
        // TODO: handle gap symbols -- replace by asterisk
        if (monomersArray[i] === gapSymbol)
          monomersArray[i] = '*';
        else
      }
      // TODO: determine conditionally (if isDna(), or isRna(), or isPeptide()) the template
      return monomersArray.join('');
    });
    return newColumn;
  }

  private convertSeparatorToFasta(): DG.Column {
    // TODO: implementation
    // * similarly to fasta2separator, divide string into monomers
    // * adjacent separators is a gap (symbol to be specified)
    // * the monomers MUST be single-character onles, otherwise forbid
    // conversion
    //getSplitterWithSeparator
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

  /** Dispatcher method for notation conversion */
  // TODO: write the bodies of converter methods
  public convert() : DG.Column {
    if (this.sourceNotation === this.targetNotation)
      throw new Error('Target notation is not specified');
    if (this.isFasta() && this.toSeparator())
      return this.convertFastaToSeparator();
    else if (this.isFasta() && this.toHelm())
      return this.convertFastaToHelm();
    else if (this.isSeparator() && this.toFasta())
      return this.convertSeparatorToFasta();
    else if (this.isSeparator() && this.toHelm())
      return this.convertSeparatorToHelm();
    else if (this.isHelm() && this.toFasta())
      return this.convertHelmToFasta();
    else
      return this.convertHelmToSeparator();
  }

  public constructor(col: DG.Column) {
    this._sourceColumn = col;
    this._sourceUnits = this._sourceColumn.tags[DG.TAGS.UNITS];
    this._sourceNotation = this.determineSourceNotation();
    this._targetNotation = this.sourceNotation;
  }
}
