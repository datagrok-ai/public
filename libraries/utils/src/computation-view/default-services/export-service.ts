/* eslint-disable valid-jsdoc */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import {ExportService, StateService} from '../common/service-interfaces';
import {tokens} from '../common/inject-tokens';

const getSheetName = (name: string, direction: DIRECTION) => {
  const idealName = `${direction} - ${name}`;
  return (idealName.length > 31) ? name.substring(0, 32) : idealName;
};

enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
}

const scalarsToSheet = (sheet: ExcelJS.Worksheet, scalars: {caption: string, value: string, units: string}[]) => {
  sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
  scalars.forEach((scalar) => {
    sheet.addRow([scalar.caption, scalar.value, scalar.units]);
  });

  sheet.getColumn(1).width = Math.max(
    ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length,
  ) * 1.2;
  sheet.getColumn(2).width = Math.max(...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
  sheet.getColumn(3).width = Math.max(...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
};

const dfToSheet = (sheet: ExcelJS.Worksheet, df: DG.DataFrame) => {
  sheet.addRow((df.columns as DG.ColumnList).names()).font = {bold: true};
  for (let i = 0; i < df.rowCount; i++)
    sheet.addRow([...df.row(i).cells].map((cell: DG.Cell) => cell.value));


  for (let i = 0; i < df.columns.length; i++) {
    sheet.getColumn(i+1).width =
      Math.max(
        ...(df.columns as DG.ColumnList).byIndex(i).categories.map((category) => category.toString().length),
        (df.columns as DG.ColumnList).byIndex(i).name.length,
      ) * 1.2;
  }
};

export class DefaultExportService implements ExportService {
  public static inject = [
    tokens.stateService,
  ] as const;

  constructor(
    private state: StateService,
  ) {}

  filename(format: string): string {
    return `${this.state.name} - ${new Date().toLocaleString()}.${this.supportedExportExtensions[format]}`;
  }

  /** Override to provide supported export formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportFormats(): string[] {
    return [
      'Excel',
    ];
  }

  /** Override to provide custom file extensions for exported formats.
   * These formats are available under the "Export" popup on the ribbon panel. */
  get supportedExportExtensions(): Record<string, string> {
    return {
      'Excel': 'xlsx',
    };
  }

  /** Override to provide custom export. */
  async export(format: string): Promise<Blob> {
    if (!this.state.func) throw new Error(`No function is associated with view`);

    const lastCall = this.state.lastCall;
    if (!lastCall) throw new Error(`Function was not called`);

    if (format !== 'Excel') throw new Error(`Format "${format}" is not supported.`);

    const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
    const exportWorkbook = new ExcelJS.Workbook();

    const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

    const isDataFrame = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

    const dfInputs = this.state.func.inputs.filter((input) => isDataFrame(input.propertyType));
    const scalarInputs = this.state.func.inputs.filter((input) => isScalarType(input.propertyType));
    const dfOutputs = this.state.func.outputs.filter((output) => isDataFrame(output.propertyType));
    const scalarOutputs = this.state.func.outputs.filter((output) => isScalarType(output.propertyType));

    dfInputs.forEach((dfInput) => {
      const visibleTitle = dfInput.options.caption || dfInput.name;
      const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT));

      const currentDf = (lastCall.inputs[dfInput.name] as DG.DataFrame);
      dfToSheet(currentDfSheet, currentDf);
    });

    if (scalarInputs.length) {
      const inputScalarsSheet = exportWorkbook.addWorksheet('Input scalars');
      scalarsToSheet(inputScalarsSheet, scalarInputs.map((scalarInput) => ({
        caption: scalarInput.options['caption'] || scalarInput.name,
        value: lastCall.inputs[scalarInput.name],
        units: scalarInput.options['units'] || '',
      })));
    }

    dfOutputs.forEach((dfOutput) => {
      const visibleTitle = dfOutput.options.caption || dfOutput.name;
      const currentDfSheet = exportWorkbook.addWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT));

      const currentDf = (lastCall.outputs[dfOutput.name] as DG.DataFrame);
      dfToSheet(currentDfSheet, currentDf);
    });

    if (scalarOutputs.length) {
      const outputScalarsSheet = exportWorkbook.addWorksheet('Output scalars');
      scalarsToSheet(outputScalarsSheet, scalarOutputs.map((scalarOutput) => ({
        caption: scalarOutput.options['caption'] || scalarOutput.name,
        value: lastCall.outputs[scalarOutput.name],
        units: scalarOutput.options['units'] || '',
      })));
    }

    const buffer = await exportWorkbook.xlsx.writeBuffer();

    return new Blob([buffer], {type: BLOB_TYPE});
  }
}
