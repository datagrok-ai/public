import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import wu from 'wu';
import {assure, getExcelColumnLetter} from "./utils";


/** Represents a well in the experimental plate */
export interface PlateWell {
  row: number;
  col: number;
  volume?: number;
  concentration?: number;
  [key: string]: any;
}


interface IPlateCsvImportOptions {
  rowColumn?: string;
  colColumn?: string;
  positionColumn?: string;
  volumeColumn?: string;
  concentrationColumn?: string;
  field?: string;   // used as a field name when a 12x8 table is passed; if not passed, defaults to 'value'
}


/** Represents experimental plate (typically 96-well assay plate) */
export class Plate {
  data: DG.DataFrame;
  rows: number = 8;
  cols: number = 12;

  constructor(rows: number, cols: number) {
    this.data = DG.DataFrame.create(rows * cols);
    this.rows = rows;
    this.cols = cols;
  }

  _idx(row: number, col: number): number {
    return row * this.cols + col;
  }

  //well(row: number, col: number): PlateWell { return new PlateWell(); }

  static fromTable(table: DG.DataFrame, field: string): Plate {
    const rows = table.rowCount;
    // remove row letters
    function containsRowLetters(c: DG.Column): boolean {
      return c.type == DG.TYPE.STRING && wu(DG.range(rows)).every((i) => getExcelColumnLetter(i) == c.get(i).toUpperCase());
    }
    const dataColumns = wu(table.columns).filter((c) => !containsRowLetters(c)).toArray();
    const cols = dataColumns.length;
    const plate = new Plate(rows, cols);

    const nonEmptyCols = dataColumns.filter((c) => c.stats.missingValueCount != rows);
    const toString = (nonEmptyCols.length > 0 && nonEmptyCols.findIndex((c) => c.type != nonEmptyCols[0].type) != -1);
    const type = toString ? DG.TYPE.STRING : dataColumns[0].type;

    const dataColumn = plate.data.columns.addNew(field, type);

    for (let colIdx = 0; colIdx < dataColumns.length; colIdx ++)
      if (dataColumns[colIdx].stats.missingValueCount != rows)
        for (let i = 0; i < rows; i++)
          dataColumn.set(plate._idx(i, colIdx), toString ? dataColumns[colIdx].get(i) : dataColumns[colIdx].getString(i), false);

    return plate;
  }

  static async fromCsvTableFile(csvPath: string, field: string, options?: DG.CsvImportOptions): Promise<Plate> {
    return this.fromTable(await grok.dapi.files.readCsv(csvPath, options), field);
  }

  static fromCsvTable(csv: string, field: string, options?: DG.CsvImportOptions): Plate {
    return this.fromTable(DG.DataFrame.fromCsv(csv), field);
  }

  static fromCsv(csv: string, options?: IPlateCsvImportOptions): Plate {
    const df = DG.DataFrame.fromCsv(csv);

    // TODO: tall format?
    // if (df.columns.contains('pos'))
    //   return new Plate(0, 0);

    return this.fromTable(df, options?.field ?? 'value');
  }

  /** Merges the attributes from {@link plates} into one plate. */
  static fromPlates(plates: Plate[]) {
    assure(plates.length > 0, 'Array is empty.');
    assure(plates.every(p => p.rows == plates[0].rows && p.cols == plates[0].cols), 'Plate dimensions differ.');

    return plates.reduce((p1, p2) => p1 ? p1.merge(p2) : p2.clone());
  }

  /** Generates [count] increasing integer numbers, starting with 0. */
  get wells(): IterableIterator<PlateWell> {
    const self = this; // Capture this context
    return (function* () {
      for (let row = 0; row < self.rows; row++) {
        for (let col = 0; col < self.cols; col++) {
          const well: PlateWell = {
            row: row,
            col: col,
          };
          for (const field of self.data.columns) {
            if (!field.isNone(self._idx(row, col))) {
              well[field.name] = field.get(self._idx(row, col));
            }
          }
          yield well;
        }
      }
    })();
  }

  /** Adds data from the other plate to this plate. Typically, you would apply plate layout */
  merge(plate: Plate) {
    for (const col of plate.data.columns) {
      this.data.columns.add(col.clone());
    }
    return this;
  }

  /** Deep cloning */
  clone(): Plate {
    const cloned = new Plate(this.rows, this.cols);
    cloned.data = this.data.clone();
    return cloned;
  }

  print() {
    for (const fieldCol of this.data.columns) {
      console.log();
      console.log(fieldCol.name);
      console.log(wu(DG.range(this.cols + 1)).map((col) => col == 0 ? '' : `${col}`).toArray().join('\t'));
      for (let row = 0; row < this.rows; row++) {
        console.log(getExcelColumnLetter(row) + '\t' +
          wu(DG.range(this.cols)).map((col) => fieldCol.getString(this._idx(row, col))).toArray().join(',\t'));
      }
    }
  }
}