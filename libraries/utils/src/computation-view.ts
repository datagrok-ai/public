/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import html2canvas from 'html2canvas';
import {FunctionView} from './function-view';

/**
 * Base class for handling Compute models (see https://github.com/datagrok-ai/public/blob/master/help/compute/compute.md).
 * In most cases, the computation is a simple {@link Func}
 * Extend it in cases where a behavior or UI not supported by the {@link FunctionView} is needed.
 *
 * It provides the following functionality out-of-the-box, where each section could be customized:
 * - a structured way to represent input and output parameters: {@link parameters}
 * - generic way to generate UI for inputs, outputs, and interactivity (running the model, etc)
 *   - persisting historical results to the db (via {@link parameters})
 * - export (to Excel and PDF): {@link export}
 * - easy loading of historical runs
 * - routing
 * - entering the real, measured (as opposed to predicted) values manually
 * - notifications for changing inputs, completion of computations, etc: {@link onInputChanged}
 * */
export class ComputationView extends FunctionView {
  /** Override to create bug reporting feature */
  reportBug: (() => Promise<void>) | null = null;

  /** Filename for exported files. Override for custom names */
  exportFilename(format: string): string {
    return `${this.name} - ${new Date().toLocaleString()}.${this.supportedExportExtensions[format]}`;
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

  /** Override to provide custom export. Obligatory to override if custom {@link buildIO()} or {@link buildOutputBlock()} is overriden. */
  async export(format: string): Promise<Blob> {
    const lastCall = this.lastCall;
    if (!lastCall) throw new Error(`Function was not called`);

    if (format !== 'Excel') throw new Error(`Format "${format}" is not supported.`);

    if (!this.func) throw new Error('The correspoding function is not specified');

    const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
    const exportWorkbook = new ExcelJS.Workbook();

    const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

    const isDataFrame = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

    const dfInputs = this.func.inputs.filter((input) => isDataFrame(input.propertyType));
    const scalarInputs = this.func.inputs.filter((input) => isScalarType(input.propertyType));
    const dfOutputs = this.func.outputs.filter((output) => isDataFrame(output.propertyType));
    const scalarOutputs = this.func.outputs.filter((output) => isScalarType(output.propertyType));

    const inputParams = [...lastCall.inputParams.values()] as DG.FuncCallParam[];
    const outputParams = [...lastCall.outputParams.values()] as DG.FuncCallParam[];

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

    const tabControl = this.resultsTabControl;
    if (tabControl) {
      for (const tabLabel of this.tabsLabels) {
        tabControl.currentPane = tabControl.getPane(tabLabel);
        await new Promise((r) => setTimeout(r, 100));
        if (tabLabel === 'Input') {
          for (const inputParam of inputParams.filter((inputParam) => inputParam.property.propertyType === DG.TYPE.DATA_FRAME)) {
            const nonGridViewers = (inputParam.aux['viewers'] as DG.Viewer[]).filter((viewer) => viewer.type !== DG.VIEWER.GRID);

            const dfInput = dfInputs.find((input) => input.name === inputParam.name);
            const visibleTitle = dfInput!.options.caption || inputParam.name;
            const currentDf = (lastCall.inputs[dfInput!.name] as DG.DataFrame);

            for (const [index, viewer] of nonGridViewers.entries()) {
              if (viewer.root.parentElement?.style.display === 'none') {
                this.paramGridSwitches.get(inputParam.name)?.click();
                await new Promise((r) => setTimeout(r, 50));
              }

              await plotToSheet(
                exportWorkbook,
                exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.INPUT)),
                viewer.root,
                currentDf.columns.length + 2,
                (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0
              );
            };
          }
        } else {
          for (const outputParam of outputParams.filter((outputParam) => outputParam.property.propertyType === DG.TYPE.DATA_FRAME && outputParam.property.category === tabLabel)) {
            const nonGridViewers = (outputParam.aux['viewers'] as DG.Viewer[]).filter((viewer) => viewer.type !== DG.VIEWER.GRID);

            const dfOutput = dfOutputs.find((input) => input.name === outputParam.name);
            const visibleTitle = dfOutput!.options.caption || outputParam.name;
            const currentDf = (lastCall.outputs[dfOutput!.name] as DG.DataFrame);

            for (const [index, viewer] of nonGridViewers.entries()) {
              if (viewer.root.parentElement?.style.display === 'none') {
                this.paramGridSwitches.get(outputParam.name)?.click();
                await new Promise((r) => setTimeout(r, 50));
              }

              await plotToSheet(
                exportWorkbook,
                exportWorkbook.getWorksheet(getSheetName(visibleTitle, DIRECTION.OUTPUT)),
                viewer.root,
                currentDf.columns.length + 2,
                (index > 0) ? Math.ceil(nonGridViewers[index-1].root.clientHeight / 20) + 1 : 0
              );
            }
          };
        }
      };
    }


    const buffer = await exportWorkbook.xlsx.writeBuffer();

    return new Blob([buffer], {type: BLOB_TYPE});
  }

  buildRibbonMenu() {
    if (this.supportedExportFormats.length === 0 && !this.reportBug && !this.helpUrl) return;

    const ribbonMenu = this.ribbonMenu.group('Model');

    if (this.supportedExportFormats.length > 0) {
      ribbonMenu
        .group('Export')
        .items(this.supportedExportFormats, async (format: string) => DG.Utils.download(this.exportFilename(format), await this.export(format)))
        .endGroup();
    }

    if (this.reportBug)
      ribbonMenu.item('Report a bug', () => this.reportBug!());

    if (this.helpUrl)
      ribbonMenu.item('Help', () => window.open(this.helpUrl!));
  }

  /** Override to create a custom input control. */
  buildRibbonPanels(): HTMLElement[][] {
    const newRibbonPanels = [
      ...this.getRibbonPanels(),
      [...(this.supportedExportFormats.length > 0) ? [ui.divH([
        ui.comboPopup(
          ui.iconFA('arrow-to-bottom'),
          this.supportedExportFormats,
          async (format: string) => DG.Utils.download(this.exportFilename(format), await this.export(format))),
      ])]: []]
    ];
    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }
}

const getSheetName = (name: string, direction: DIRECTION) => {
  const idealName = `${direction} - ${name}`;
  return (idealName.length > 31) ? name.substring(0, 32) : idealName;
};

enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
}

const scalarsToSheet = (sheet: ExcelJS.Worksheet, scalars: { caption: string, value: string, units: string }[]) => {
  sheet.addRow(['Parameter', 'Value', 'Units']).font = {bold: true};
  scalars.forEach((scalar) => {
    sheet.addRow([scalar.caption, scalar.value, scalar.units]);
  });

  sheet.getColumn(1).width = Math.max(
    ...scalars.map((scalar) => scalar.caption.toString().length), 'Parameter'.length
  ) * 1.2;
  sheet.getColumn(2).width = Math.max(...scalars.map((scalar) => scalar.value.toString().length), 'Value'.length) * 1.2;
  sheet.getColumn(3).width = Math.max(...scalars.map((scalar) => scalar.units.toString().length), 'Units'.length) * 1.2;
};

const dfToSheet = (sheet: ExcelJS.Worksheet, df: DG.DataFrame) => {
  sheet.addRow((df.columns as DG.ColumnList).names()).font = {bold: true};
  for (let i = 0; i < df.rowCount; i++)
    sheet.addRow([...df.row(i).cells].map((cell: DG.Cell) => cell.value));

  for (let i = 0; i < df.columns.length; i++) {
    sheet.getColumn(i + 1).width =
      Math.max(
        ...(df.columns as DG.ColumnList).byIndex(i).categories.map((category) => category.toString().length),
        (df.columns as DG.ColumnList).byIndex(i).name.length
      ) * 1.2;
  }
};

const plotToSheet = async (exportWb: ExcelJS.Workbook, sheet: ExcelJS.Worksheet, plot: HTMLElement, columnForImage: number, rowForImage: number = 0) => {
  const canvas = await html2canvas(plot as HTMLElement, {logging: false});
  const dataUrl = canvas.toDataURL('image/png');

  const imageId = exportWb.addImage({
    base64: dataUrl,
    extension: 'png',
  });
  sheet.addImage(imageId, {
    tl: {col: columnForImage, row: rowForImage},
    ext: {width: canvas.width, height: canvas.height},
  });
};
