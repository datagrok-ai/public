/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import ExcelJS from 'exceljs';
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
  _inputFields: Map<string, DG.InputBase> = new Map<string, DG.InputBase>();

  onViewInitialized: rxjs.Subject<void> = new rxjs.Subject();
  onInputChanged: rxjs.Subject<void> = new rxjs.Subject();

  onComputationStarted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted on computation start
  onComputationCompleted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly or not
  onComputationSucceeded: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly
  onComputationError: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended with error

  constructor(func: DG.Func) {
    super(func);
  }

  /** Override if custom initialization is required */
  init() {
    super.init();
    this.renderExportPanel();
    this.renderHistoryPanel();
    this.renderRibbonPanels();

    this.onComputationSucceeded.subscribe(async () => {
      await this.renderOutputPanel();
    });

    this.onComputationError.subscribe(() => {
      grok.shell.error('Computation failed');
    });
  }

  /** All inputs that are bound to fields */
  get inputFields(): Map<string, DG.InputBase> {return this._inputFields;}

  /** Panel for exporting functions file */
  exportPanel: HTMLElement = ui.div();
  /** Override if custom export panel is required */
  buildExportPanel() {
    return ui.comboPopup(
      ui.icons.save(null, 'Export'),
      this.supportedExportFormats,
      async (format: string) =>
        DG.Utils.download(this.exportFilename(format), await this.export(format)),
    );
  }
  /** Call to update export panel */
  renderExportPanel() {
    const newExportPanel = this.buildExportPanel();
    this.exportPanel.replaceWith(newExportPanel);
    this.exportPanel = newExportPanel;
  }

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

  /** Override to provide custom export. */
  async export(format: string): Promise<Blob> {
    if (!this.func) throw new Error(`No function is associated with view`);

    const lastCall = this.lastCall;
    if (!lastCall) throw new Error(`Function was not called`);

    if (format !== 'Excel') throw new Error(`Format "${format}" is not supported.`);

    const BLOB_TYPE = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;charset=UTF-8';
    const exportWorkbook = new ExcelJS.Workbook();

    const isScalarType = (type: DG.TYPE) => (DG.TYPES_SCALAR.has(type));

    const isDataFrame = (type: DG.TYPE) => (type === DG.TYPE.DATA_FRAME);

    const dfInputs = this.func.inputs.filter((input) => isDataFrame(input.propertyType));
    const scalarInputs = this.func.inputs.filter((input) => isScalarType(input.propertyType));
    const dfOutputs = this.func.outputs.filter((output) => isDataFrame(output.propertyType));
    const scalarOutputs = this.func.outputs.filter((output) => isScalarType(output.propertyType));

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

  /** Panel to deal with historical runs */
  historyPanel: HTMLElement = ui.div();
  /** Override to provide custom historical runs UI */
  buildHistoryPanel() {
    return ui.iconFA('history', () => {
      this.pullRuns().then(async (historicalRuns) => {
        const menu = DG.Menu.popup();
        let i = 0;
        for (const run of historicalRuns) {
          //@ts-ignore
          menu.item(`${run.func.friendlyName} — ${i}`, async () => this.historicalRunService.loadRun(run.id));
          i++;
        }
        menu.show();
      });
    });
  }
  /** Call to update history panel */
  renderHistoryPanel() {
    const newHistoryPanel = this.buildHistoryPanel();
    this.historyPanel.replaceWith(newHistoryPanel);
    this.historyPanel = newHistoryPanel;
  }

  /** Saves the computation results to the historical results, returns its id. See also {@link loadRun}. */
  async saveRun(call: DG.FuncCall): Promise<string> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.save(call);*/
  }
  /** Loads the specified historical results. See also {@link saveRun}. */
  async loadRun(runId: string): Promise<DG.FuncCall> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.find(call.id);*/
  }
  /** Loads all the function call of this function. */
  async pullRuns(): Promise<DG.FuncCall[]> {
    throw new Error('Method is not implemented');
    /* await grok.dapi.functions.calls.list();*/
  }

  /** Pass function arguments' values before computing */
  async prepare() {
    if (!this.func) return;

    this.lastCall = this.func.prepare()!;
  }

  /** Maps inputs to parameters, computes, and maps output parameters to the UI. */
  async run(): Promise<void> {
    await this.prepare();
    if (!this.lastCall) throw new Error('Func is not prepared!');

    // this.inputFieldsToParameters(this.lastCall);

    try {
      await this.compute(this.lastCall);
    } catch (e) {
      this.onComputationError.next(this.lastCall);
    }

    console.log('complteed', this.lastCall);
    this.onComputationCompleted.next(this.lastCall);

    // this.outputParametersToView(this.lastCall);
  }

  /** Override to customize ribbon panels */
  public async buildRibbonPanels(): Promise<HTMLElement[][]> {
    return [[
      this.buildExportPanel(),
    ]];
  }
  /** Call to update ribbon panels */
  public async renderRibbonPanels() {
    this.setRibbonPanels(await this.buildRibbonPanels());
  }

  /** Builds the complete view. Override it if completely custom implementation is needed,
   * otherwise consider overriding {@link buildInputPanel} and/or {@link buildOutputPanel} */
  override build(): HTMLElement {
    return super.build();
  }

  /** Override to build completely custom input block (the one with inputs, usually on the left).
   * If only specific field, consider overriding {@link buildCustomInputs}. */
  override buildInputPanel(): HTMLElement {
    return super.buildInputPanel();
  }

  /** Custom inputs for the specified fields. Inputs for these fields will not be created by {@link buildInputPanel}. */
  buildCustomInputs(): Map<String, DG.InputBase> {return new Map();}

  /** Override to create output block. */
  override buildOutputPanel(): HTMLElement {
    console.log('ouptut updated');
    const outputs = this.lastCall?.outputs as Record<string, any>;
    const panel = ui.accordion('Output data');
    panel.addPane('Output data', () => {
      return ui.divV(
        Object.entries(outputs).map(([key, val]) => ui.span([`${key}: `, `${val}`])),
      );
    });
    return panel.root;
  }

  /** Override to create a custom input control. */
  buildRibbonPanel(): HTMLElement {
    return ui.divH([
      ui.button('RUN', () => {}),
      ui.comboPopup('Export',
        this.supportedExportFormats, (format: string) => this.export(format)),
    ]);
  }

  /** Creates function parameters based on the data entered by the user. */
  inputFieldsToParameters(call: DG.FuncCall): void {super.inputFieldsToParameters(call);}

  /** Visualizes computation results */
  outputParametersToView(call: DG.FuncCall): void {super.outputParametersToView(call);}

  /** Override to provide custom computation error handling. */
  processComputationError(call: DG.FuncCall) {grok.shell.error(call.errorMessage!);}

  /** helper methods */
  override buildInputForm(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        await this.run();
      });
      const editor = ui.div([], 'ui-form');
      const inputs: DG.InputBase[] = await call.buildEditor(editor, {condensed: true});
      editor.appendChild(ui.divH([ui.divV([this.historyPanel], {style: {'justify-content': 'center'}}), runButton], {style: {'justify-content': 'space-between'}}));
      return ui.panel([editor]);
    });
  };
}

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
