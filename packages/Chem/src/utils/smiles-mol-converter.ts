/* Do not change these import lines to match external modules in webpack configuration */
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// import {getRdKitModule} from './chem-common-rdkit';

/** Class for handling convresion from SMILES to Molfile and vice versa */
export class NotationConverter {
  private readonly _column: DG.Column; // the column to be converted
  private readonly _units: string; // units, smiles/molblock

  /**
   * Create a new empty column of the specified notation type and the same
   * length as this._column
   *
   * @param {string} targetNotation
   * @return {DG.Column}
   */
  private getNewColumn(targetNotation: string): DG.Column {
    const col = this._column;
    const len = col.length;
    const name = targetNotation.toLowerCase() + '(' + col.name + ')';
    const newColName = col.dataFrame.columns.getUnusedName(name);
    const newColumn = DG.Column.fromList('string', newColName, new Array(len).fill(''));
    newColumn.semType = DG.SEMTYPE.MOLECULE;
    newColumn.setTag(DG.TAGS.UNITS, targetNotation);
    newColumn.setTag(DG.TAGS.CELL_RENDERER, 'Molecule');

    return newColumn;
  }
  public async convertSmilesToV2000(): Promise<DG.Column> {
    // const mod = await getRdKitModule();
    const targetNotation = DG.UNITS.Molecule.MOLBLOCK;

    // assign the values to the empty column
    const newColumn = this.getNewColumn(targetNotation);
    newColumn.init((idx: number) => {
      // const smilesString = this._column.get(idx);
      return 'dummy';
    });
    return newColumn;
  };

  // public convertSmilesToV3000(): DG.Column {};
  // public convertV2000ToV3000(): DG.Column {};
  // public convertV2000ToSmiles(): DG.Column {};
  // public convertV3000ToV2000(): DG.Column {};
  // public convertV3000ToSmiles(): DG.Column {};

  public constructor(col: DG.Column) {
    this._column = col;
    const units = this._column.tags[DG.TAGS.UNITS];
    if (units == DG.UNITS.Molecule.MOLBLOCK || units == DG.UNITS.Molecule.SMILES)
      this._units = units;
    else
      throw new Error(`Units specified in the colum are: ${units}`);
  }
}
