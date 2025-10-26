import * as DG from 'datagrok-api/dg';
import * as ExcelJS from 'exceljs';
import wu from 'wu';
import {Plate} from './plate';
import {DataFrame} from 'datagrok-api/dg';

interface ExcelPlatePosition {
  sheet: ExcelJS.Worksheet;
  startRow: number;
  startCol: number;
  rows: number;
  cols: number;
}


/** Finds plate positions in the Excel workbook by iterating over sheets
 * and finding row/column label patterns such as (1-24, A-P). */
export function findPlatePositions(workbook: ExcelJS.Workbook): ExcelPlatePosition[] {
  const results: ExcelPlatePosition[] = [];
  const ACharCode = 'A'.charCodeAt(0);
  workbook.eachSheet((worksheet) => {
    // TODO: Hack for current things, remove later
    if (worksheet.name?.toLowerCase() === 'desired outputs')
      return;
    // Create a grid representation of the worksheet
    const grid: any[][] = [];
    worksheet.eachRow((row, rowIndex) => {
      grid[rowIndex] = [];
      row.eachCell((cell, colIndex) => {
        grid[rowIndex][colIndex] = cell.value;
      });
    });

    // Find numeric column headers (1, 2, 3, ...)
    let headerRow = -1; let startCol = -1; let numColumns = 0;
    for (let i = 1; i < grid.length; i++) {
      const row = grid[i] || [];
      let numSequence = 0; let seqStart = -1;

      for (let j = 1; j < row.length; j++) {
        if (typeof row[j] === 'number' && Number.isInteger(row[j]) && row[j] === numSequence + 1) {
          if (numSequence === 0) seqStart = j;
          numSequence++;
        } else if (numSequence > 0) { break; }
      }

      if (numSequence >= 12) { // Minimum 12 columns
        headerRow = i;
        startCol = seqStart;
        numColumns = numSequence;
        break;
      }
    }

    // Find alphabetic row headers (A, B, C, ...)
    let headerCol = -1; let startRow = -1; let numRows = 0;
    if (headerRow !== -1) {
      for (let j = 1; j < startCol; j++) {
        let letterSequence = 0; let seqStart = -1;

        for (let i = headerRow + 1; i < grid.length; i++) {
          if (!grid[i]) continue;
          const val = grid[i][j];

          if (typeof val === 'string' && /^[A-Za-z]{1,2}$/.test(val)) {
            const prefix = letterSequence > 25 ? 'A' : '';
            const expected = prefix + String.fromCharCode(ACharCode + letterSequence % 26);
            if (val?.toLowerCase() === expected?.toLowerCase()) {
              if (letterSequence === 0) seqStart = i;
              letterSequence++;
            } else { break; }
          } else if (letterSequence > 0) { break; }
        }

        if (letterSequence >= 8) { // Minimum 8 rows
          headerCol = j;
          startRow = seqStart;
          numRows = letterSequence;
          break;
        }
      }
    }

    if (headerRow !== -1 && headerCol !== -1) {
      results.push({
        sheet: worksheet,
        startRow: headerRow,
        startCol: headerCol,
        rows: numRows,
        cols: numColumns
      });
    }
  });

  return results;
}


/** Constructs a plate from the specified position in the Excel sheet. */
/** Constructs a plate from the specified position in the Excel sheet. */
export function getPlateFromSheet(p: ExcelPlatePosition): Plate {
  const df = DataFrame.create(p.rows);

  for (let col = p.startCol + 1; col <= p.startCol + p.cols; col++) {
    const values: any[] = [];
    let hasNonNumeric = false;

    // First pass: collect values and check if they're all numeric
    for (let i = 0; i < p.rows; i++) {
      const cellValue = p.sheet.getCell(p.startRow + i + 1, col).value;
      values.push(cellValue);

      if (cellValue !== null && cellValue !== undefined && cellValue !== '' &&
          typeof cellValue !== 'number' && isNaN(Number(cellValue)))
        hasNonNumeric = true;
    }

    // Create appropriate column type
    const colName = `${col - p.startCol}`;
    if (hasNonNumeric) {
      // String column
      const stringCol = DG.Column.fromStrings(colName, values.map((v) => v?.toString() ?? ''));
      df.columns.add(stringCol);
    } else {
      // Numeric column
      const numericCol = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, colName,
        values.map((v) => {
          if (v === null || v === undefined || v === '') return null;
          return typeof v === 'number' ? v : Number(v);
        })
      );
      df.columns.add(numericCol);
    }
  }

  return Plate.fromGridDataFrame(df, p.sheet.name);
}
