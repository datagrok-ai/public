/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SeqHandler} from './seq-handler';
import {NOTATION} from './macromolecule/consts';

/** Class for parsing FASTA files */
export class FastaFileHandler {
  private _fileContent: string;
  private _descriptionsArray: string[] = []; // parsed FASTA descriptions
  private _sequencesArray: string[] = []; // parsed FASTA sequeces
  // private _columnsParsed: boolean = false;

  public get descriptionsArray(): string[] { return this._descriptionsArray; }

  public get sequencesArray(): string[] { return this._sequencesArray; }

  /**
   * Helper method to parse a macromolecule from a FASTA file (string)
   *
   * @param {number} startOfSequence  index of macromolecule substring beginning
   * @param {number} endOfSequence  index of macromolecule substring end

   * @return {string} parsed macromolecule
   */
  private parseMacromolecule(
    startOfSequence: number,
    endOfSequence: number
  ): string {
    const seq = this._fileContent.slice(startOfSequence, endOfSequence);
    const seqArray = seq.split(/\s/);
    return seqArray.join('');
  }

  /** Parse descriptions and sequences from a FASTA string */
  private parseColumns() {
    const regex = /^>(.*)$/gm; // match 'description' lines starting with >

    let startOfSequence = 0;
    let match; // match.index is the beginning of the matched line
    while (match = regex.exec(this._fileContent)) {
      const description = this._fileContent.substring(match.index + 1, regex.lastIndex);
      this._descriptionsArray.push(description);
      if (startOfSequence !== 0)
        this._sequencesArray.push(this.parseMacromolecule(startOfSequence, match.index));
      startOfSequence = regex.lastIndex + 1;
    }
    this._sequencesArray.push(this.parseMacromolecule(startOfSequence, -1));

    // this._columnsParsed = true;
  }

  /**
   * File-handler method for import as FASTA
   *
   * @return {DG.DataFrame[]} dataframe with parsed FASTA content
   */
  public importFasta(): DG.DataFrame [] {
    const descriptionsArrayCol = DG.Column.fromStrings('description', this.descriptionsArray);
    const sequenceCol = DG.Column.fromStrings('sequence', this.sequencesArray);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.meta.units = NOTATION.FASTA;

    // here should go the code from units handler
    const sh = SeqHandler.forColumn(sequenceCol);

    return [DG.DataFrame.fromColumns([
      descriptionsArrayCol,
      sequenceCol,
    ])];
  }

  constructor(fileContent: string) {
    this._fileContent = fileContent;
    this.parseColumns();
  }
}
