/* eslint-disable max-len */
/* eslint-disable valid-jsdoc */
import wu from 'wu';
import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';

export class FunctionView extends DG.ViewBase {
  onViewInitialized: rxjs.Subject<void> = new rxjs.Subject();
  onInputChanged: rxjs.Subject<void> = new rxjs.Subject();

  onComputationStarted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted on computation start
  onComputationCompleted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly or not
  onComputationSucceeded: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly
  onComputationError: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended with error

  constructor(func: DG.Func) {
    super();
    this.func = func;
  }

  /** Override if custom initialization is required */
  init() {
    this.lastCall ??= this.func.prepare();
    this.lastCall.aux['view'] = this;
    this.lastCall.context = DG.Context.cloneDefault();

    this.singleDfParam = wu(this.func.outputs).filter((p) => p.propertyType == DG.TYPE.DATA_FRAME)
      .toArray().length == 1;

    for (const inParam of wu(this.lastCall.inputParams.values() as DG.FuncCallParam[])
      .filter((p: DG.FuncCallParam) =>
        p.property.propertyType == DG.TYPE.DATA_FRAME && p.property.options['viewer'] != null)) {
      this.showInputs = true;
      const self = this;
      this.subs.push(inParam.onChanged.subscribe(async function(param: DG.FuncCallParam) {
        self.clearResults(true, false);
        param.processOutput();

        self.appendResultDataFrame(param,
          {height: (self.singleDfParam && !grok.shell.windows.presentationMode) ? 600 : 400,
            category: 'INPUT'});
      }));

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

    this.root.appendChild(this.build());
    this.name = this.func.friendlyName;
    this.box = true;
  }

  private _type: string = 'function';
  public get type(): string {
    return this._type;
  }
  public func: DG.Func;
  public lastCall: DG.FuncCall | null = null;
  private _inputFields: Map<string, DG.InputBase> = new Map<string, DG.InputBase>();

  /** All inputs that are bound to fields */
  get inputFields(): Map<string, DG.InputBase> {return this._inputFields;}

  /** Input panel for form, run button, etc. On the left side by default */
  inputPanel = ui.div() as HTMLElement;

  /** Override to build completely custom input block (the one with inputs, usually on the left).
   * If only specific field, consider overriding {@link buildCustomInputs}. */
  buildInputPanel(): HTMLElement {
    return ui.div([this.buildInputForm(this.lastCall!)]);
  }

  /** Call to update input panel */
  async renderInputPanel() {
    const newInputPanel = this.buildInputPanel();
    this.inputPanel.replaceWith(newInputPanel);
    this.inputPanel = newInputPanel;
  }

  /** Output panel for viewing results. On the right side by default */
  outputPanel = ui.div() as HTMLElement;
  /** Override if custom ouptut panel is required */
  buildOutputPanel(): HTMLElement {
    return ui.div();
  }
  /** Call to update output panel */
  async renderOutputPanel() {
    const newOutputPanel = this.buildOutputPanel();
    console.log('newOutputPanel', newOutputPanel);
    this.outputPanel.replaceWith(newOutputPanel);
    this.outputPanel = newOutputPanel;
    console.log('outputPanel', this.outputPanel);
  }

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
  build(): HTMLElement {
    this.renderInputPanel();
    this.renderOutputPanel();
    return ui.divH([this.inputPanel, this.outputPanel]);
  }

  public async render() {
    const newRoot = this.build();
    this.root.replaceWith(newRoot);
    this.root = newRoot;
  }

  /** Override to create a custom input control. */
  buildRibbonPanel(): HTMLElement {
    return ui.divH([
      ui.button('RUN', () => {}),
      ui.comboPopup('Export',
        this.supportedExportFormats, (format: string) => this.export(format)),
    ]);
  }

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

  async compute(call: DG.FuncCall): Promise<void> {
    ui.setUpdateIndicator(this.outputPanel, true);
    this.lastCall = await call.call(true, undefined, {processed: true});
    ui.setUpdateIndicator(this.outputPanel, false);
  }

  showInputs: boolean = false;
  singleDfParam: boolean = false;
  paramViewers: Map<string, DG.Viewer[]> = new Map();
  resultsTabControl: DG.TabControl | undefined;
  resultTabs: Map<String, HTMLElement> = new Map();
  resultsDiv: HTMLElement = ui.panel([], 'grok-func-results');
  inputsDiv: HTMLDivElement = ui.panel([], 'grok-func-results');

  inputFieldsToParameters(call: DG.FuncCall): void { }

  /** helper methods */
  appendResultDataFrame(param: DG.FuncCallParam, options?: { caption?: string, category?: string, height?: number}) {
    const df = param.value;
    let caption = options?.caption;
    let height = options?.height ?? 400;
    const viewers: DG.Viewer[] = param.aux['viewers'] ?? [];
    caption ??= param.aux['viewerTitle'] ??
      ((this.singleDfParam && viewers.length == 1) ? '' : param.property.caption) ?? '';
    let existingViewers: DG.Viewer[] | undefined = this.paramViewers.get(param.name);
    if (existingViewers != null) {
      for (const v of existingViewers)
        v.dataFrame = df;

      return;
    }
    existingViewers ??= [];
    this.paramViewers.set(param.name, existingViewers);
    if (viewers.length == 0)
      viewers.push(df.plot.grid());

    const hideList: HTMLElement[] = [];
    let blocks: number = 0;
    let blockSize: number = 0;
    const gridWrapper = ui.block([]);

    const gridSwitch = !wu(viewers).some((v: DG.Viewer) => v.type == 'Grid');
    $(gridWrapper).hide();

    const getHeader = (sw: boolean) => {
      const s = caption ?? df.name ?? '';
      const header = ui.div([], 'grok-func-results-header');
      if (s != '') {
        const h = ui.h1(s);
        ui.Tooltip.bind(h, () => s);
        header.appendChild(h);
      }
      if (gridSwitch) {
        const icon = ui.iconSvg('table', (e) => {
          e.stopPropagation();
          ui.setDisplayAll(hideList, !sw);
          ui.setDisplay(gridWrapper, sw);
          gridWrapper.classList.add(`ui-block-${blockSize}`);
        }
        , 'Show grid');
        if (!sw)
          icon.classList.add('active');
        header.appendChild(icon);
      }
      if (!grok.shell.tables.includes(df)) {
        header.appendChild(ui.icons.add((e: any) => {
          e.stopPropagation();
          const v = grok.shell.addTableView(df);
          (async () => {
            for (const viewer of viewers) {
              if (viewer.type != 'Grid') {
                const newViewer = await df.plot.fromType(viewer.type) as DG.Viewer;
                newViewer.setOptions(viewer.getOptions());
                v.addViewer(newViewer);
              }
            }
          })();
        },
        'Add to workspace',
        ));
      }
      return header;
    };

    if (gridSwitch) {
      gridWrapper.appendChild(getHeader(false));
      const grid = df.plot.grid();
      grid.root.style.height = `${height}px`;
      gridWrapper.appendChild(grid.root);
      existingViewers.push(grid);
    }

    let header = getHeader(true);
    this._appendResultElement(gridWrapper, options?.category);
    for (const viewer of viewers) {
      if (viewer?.tags == null)
        continue;
      existingViewers.push(viewer);
      const block = viewer.tags['.block-size'] ?? 100;
      const wrapper = ui.block([], `ui-block-${block}`);
      if (blocks + block <= 100) {
        hideList.push(wrapper);
        blockSize += block;
      }
      blocks += block;
      wrapper.appendChild(header);
      header = ui.div([], 'grok-func-results-header');
      wrapper.appendChild(viewer.root);
      if (viewer.type == 'grid') {
        // @ts-ignore
        const totalHeight = viewer.getOptions()['rowHeight'] *
          (viewer.dataFrame!.rowCount + 1);
        if (totalHeight < height)
          height = totalHeight;
      }
      viewer.root.style.height = `${height}px`;
      this._appendResultElement(wrapper, options?.category);
    }
  }

  _appendResultElement(d: HTMLElement, category?: string) {
    if (category != null && this.resultsTabControl != undefined && this.resultTabs.get(category) != null)
      this.resultTabs.get(category)!.appendChild(d);
    else
      this.resultsDiv.appendChild(d);
  }

  /** helper methods */
  buildInputForm(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        call.aux['view'] = this.dart;
        await this.run();
      });
      const editor = ui.div();
      const inputs: DG.InputBase[] = await call.buildEditor(editor, {condensed: true});
      editor.appendChild(ui.buttonsInput([runButton]));
      for (const input of inputs)
        this._inputFields.set(input.property.name, input);
      return editor;
    });
  };

  clearResults(clearTabs: boolean = true, showOutput: boolean = true) {
    const categories: string[] = [];
    //  resultTabs.clear();
    if (this.showInputs) {
      categories.push('INPUT');
      this.resultTabs.set('INPUT', this.inputsDiv);
    }
    for (const p of this.func!.outputs) {
      if (categories.includes(p.category))
        continue;
      categories.push(p.category);
    }
    if (clearTabs && (categories.length > 1 || (categories.length == 1 && categories[0] != 'Misc'))) {
      this.resultsTabControl = DG.TabControl.create();
      for (const c of categories) {
        if (!this.resultTabs.has(c))
          this.resultTabs.set(c, ui.div([], 'ui-panel,grok-func-results'));
        let name = c;
        if (this.showInputs && categories.length == 2 && c == 'Misc')
          name = 'OUTPUT';
        this.resultsTabControl.addPane(name, () => this.resultTabs.get(c) ?? ui.div());
      }
      if (categories.length > 1 && this.showInputs && showOutput)
        this.resultsTabControl.currentPane = this.resultsTabControl.panes[1];
      this.resultsDiv = this.resultsTabControl.root;
    } else {
      this.resultsDiv = ui.panel([], 'grok-func-results');
      this.resultsTabControl = undefined;
    }
    this.buildOutputPanel();
  }

  outputParametersToView(call: DG.FuncCall): void {
    this.clearResults(false, true);
    for (const [, p] of call.outputParams) {
      p.processOutput();
      if (p.property.propertyType == DG.TYPE.DATA_FRAME && p.value != null)
        this.appendResultDataFrame(p, {caption: p.property.name, category: p.property.category});
    }
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

