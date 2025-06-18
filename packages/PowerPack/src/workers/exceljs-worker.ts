import ExcelJS from 'exceljs';

onmessage = async (event: {data: {data: Uint8Array, sheetName?: string}}) => {
  const bytes = event.data.data;
  const sheetName = event.data.sheetName;
  const workbook = new ExcelJS.Workbook();
  const wb = await workbook.xlsx.load(bytes);

  const convertWorkbookToDataFrames = (workbook: ExcelJS.Workbook) => {
    const dfNames: string[] = [];
    const dfs = workbook.worksheets.filter((sheet) => !sheetName || sheet.name === sheetName).map((sheet) => {
      const rows: (string | null)[][] = convertSheetToDataFrame(sheet);
      dfNames.push(sheet.name);
      return rows;
    });
    return {result: dfs, names: dfNames};
  };

  const convertSheetToDataFrame = (sheet: ExcelJS.Worksheet): (string | null)[][] => {
    const rows: (string | null)[][] = new Array(sheet.rowCount).fill(null).map((() => []));
    sheet.eachRow({includeEmpty: true}, (row, rowIndex) => {
      const values = Array.isArray(row.values) ? row.values.slice(1) : []; // 1-based row.values, 0 index is always empty
      rows[rowIndex - 1] = values.map((cell) => (cell != null ? cell.toString() : null));
    });
    return rows;
  };

  postMessage(convertWorkbookToDataFrames(wb));
};
