/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import wu from 'wu';
import $ from 'cash-dom';
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import ExcelJS from 'exceljs';
import html2canvas from 'html2canvas';

export class FunctionView extends DG.ViewBase {
  constructor(funcCall: DG.FuncCall | null) {
    super();
    this.preparedCall = funcCall;
    this.context = DG.Context.cloneDefault();
    this.controlsRoot.style.maxWidth = '370px';
    this.box = true;

    if (funcCall)
      this.init(funcCall);
  }

  get func() {
    return this.preparedCall!.func;
  }

  onViewInitialized: rxjs.Subject<void> = new rxjs.Subject();
  onInputChanged: rxjs.Subject<void> = new rxjs.Subject();

  onComputationStarted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted on computation start
  onComputationCompleted: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly or not
  onComputationSucceeded: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended properly
  onComputationError: rxjs.Subject<DG.FuncCall> = new rxjs.Subject(); // emitted anyway if computation ended with error

  private _type: string = 'function';
  public get type(): string {
    return this._type;
  }

  public get isInputPanelRequired() {
    return this.func.inputs.some((p) => p.propertyType == DG.TYPE.DATA_FRAME && p.options['viewer'] != null);
  }

  public get tabsLabels() {
    return Object.keys(this.paramToCategoryMap);
  }

  readonly context: DG.Context;
  public preparedCall: DG.FuncCall | null;
  public lastCall: DG.FuncCall | null = null;
  protected outputCategories = [] as string[];
  protected paramToCategoryMap = {} as Record<string, string[]>;

  async init(funcCall: DG.FuncCall) {
    this.preparedCall ??= funcCall;

    this.preparedCall.aux['view'] = this;
    this.preparedCall.context = this.context;

    this.singleDfParam = wu(this.func.outputs).filter((p) => p.propertyType == DG.TYPE.DATA_FRAME)
      .toArray().length == 1;

    for (const inParam of wu(this.preparedCall.inputParams.values() as DG.FuncCallParam[])
      .filter((p: DG.FuncCallParam) =>
        p.property.propertyType == DG.TYPE.DATA_FRAME && p.property.options['viewer'] != null)) {
      this.showInputs = true;
      const self = this;
      this.subs.push(inParam.onChanged.subscribe(async function(param: DG.FuncCallParam) {
        self.clearResults(false);
        param.processOutput();

        self.appendResultDataFrame(param, {
          height: (self.singleDfParam && !grok.shell.windows.presentationMode) ? 600 : 400,
          category: 'Input'
        });
      }));
    }

    const outParamCategories = [
      ...new Set(
        this.func!.outputs
          .map((p) => p.category) // get all output params' categories
          .filter((category) => category !== 'Misc') // filter out the default ones
      )]; // keep only unique of them

    this.outputCategories = [
      ...!outParamCategories.length ? ['Output'] : [], // if no categories are stated, the default category is added
      ...outParamCategories,
    ];

    if (this.isInputPanelRequired)
      this.func.inputs.forEach((p) => this.paramToCategoryMap['Input'] ? this.paramToCategoryMap['Input'].push(p.name): this.paramToCategoryMap['Input'] = [p.name]);

    this.func.outputs.forEach((p) => this.paramToCategoryMap[p.category === 'Misc' ? 'Output': p.category] ? this.paramToCategoryMap[p.category === 'Misc' ? 'Output': p.category].push(p.name) : this.paramToCategoryMap[p.category === 'Misc' ? 'Output': p.category] = [p.name]);

    this.root.appendChild(this.build());

    this.buildRibbonPanels();
    this.buildHistoryBlock();
    this.buildRibbonMenu();

    this.name = this.func.friendlyName;
  }

  /** Override to create a fully custom UI */
  build(): HTMLElement {
    const inputBlock = this.buildInputBlock();
    const outputBlock = this.buildOutputBlock();

    return ui.splitH([inputBlock, outputBlock]);
  }

  /** Override to create a custom input block. */
  buildInputBlock(): HTMLElement {
    const funcDiv = ui.div([this.renderRunSection(this.preparedCall!)], 'ui-div');
    this.controlsRoot.innerHTML = '';
    this.controlsRoot.appendChild(funcDiv);
    return this.controlsRoot;
  }

  /** Override to create a custom output block. */
  buildOutputBlock(): HTMLElement {
    this.resultsRoot.innerHTML = '';
    this.resultsRoot.appendChild(this.resultsDiv);
    return this.resultsRoot;
  }

  /** Override to create a custom historical runs control. */
  buildHistoryBlock(): HTMLElement {
    const newHistoryBlock = ui.iconFA('history', () => {
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
    this.historyRoot.innerHTML = '';
    this.historyRoot.append(newHistoryBlock);
    return newHistoryBlock;
  }

  /** Override to create a custom input control. */
  buildRibbonPanels(): HTMLElement[][] {
    const newRibbonPanels = [...this.getRibbonPanels(), [ui.divH([
      ui.comboPopup(
        ui.iconFA('arrow-to-bottom'),
        this.supportedExportFormats,
        async (format: string) => DG.Utils.download(this.exportFilename(format), await this.export(format))),
    ])]];
    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }

  /** Override to create a custom ribbon menu on the top. */
  buildRibbonMenu() {
    const ribbonMenu = this.ribbonMenu
      .group('Model')
      .group('Export')
      .items(this.supportedExportFormats, async (format: string) => DG.Utils.download(this.exportFilename(format), await this.export(format)))
      .endGroup();

    if (this.reportBug)
      ribbonMenu.item('Report a bug', () => this.reportBug!());

    if (this.helpUrl)
      ribbonMenu.item('Get help', () => window.open(this.helpUrl!));
  }

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

  /** Override to provide custom export. Obligatory to override if custom {@link build()} or {@link buildOutputBlock()} is overriden. */
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


  clearResults(switchToOutput: boolean = true) {
    const categories: string[] = [];
    // this.resultTabs.clear();
    if (this.showInputs) {
      categories.push('Input');
      this.resultTabs.set('Input', this.inputsDiv);
    }
    for (const p of this.func.outputs) {
      if (categories.includes(p.category))
        continue;
      categories.push(p.category);
    }
    if ((categories.length > 1 || (categories.length == 1 && categories[0] != 'Misc'))) {
      this.resultsTabControl = DG.TabControl.create();
      for (const c of categories) {
        if (!this.resultTabs.has(c))
          this.resultTabs.set(c, ui.div([], 'ui-panel, grok-func-results, ui-box'));
        let name = c;
        if (this.showInputs && categories.length == 2 && c == 'Misc')
          name = 'OUTPUT';
        this.resultsTabControl.addPane(name, () => this.resultTabs.get(c) ?? ui.div());
      }
      if (categories.length > 1 && this.showInputs && switchToOutput)
        this.resultsTabControl.currentPane = this.resultsTabControl.panes[1];
      this.resultsDiv = this.resultsTabControl.root;
    }
  }

  async run(): Promise<void> {
    try {
      await this.compute(this.preparedCall!);
    } catch (e) {
      // this.onComputationError.next(this.lastCall);
    }
    this.onComputationSucceeded.next(this.lastCall!);
    this.lastCall = this.preparedCall;

    this.outputParametersToView(this.lastCall!);
  }

  async compute(call: DG.FuncCall): Promise<void> {
    await call.call(true, undefined, {processed: true});
  }

  outputParametersToView(call: DG.FuncCall): void {
    this.clearResults(true);
    for (const p of call.outputParams.values() as DG.FuncCallParam[]) {
      p.processOutput();
      if (p.property.propertyType == DG.TYPE.DATA_FRAME && p.value != null)
        this.appendResultDataFrame(p, {caption: p.property.name, category: p.property.category});
      else
        this.appendResultScalar(p, {caption: p.property.name, category: p.property.category});
    }
    this.buildOutputBlock();
  }

  renderRunSection(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        call.aux['view'] = this.dart;
        await this.run();
      });
      const editor = ui.div();
      const inputs: DG.InputBase[] = await call.buildEditor(editor, {condensed: true});
      editor.classList.add('ui-form');
      const buttons = ui.divH([this.historyRoot, runButton], {style: {'justify-content': 'space-between'}});
      editor.appendChild(buttons);
      return editor;
    });
  }

  // flag to show or not the "Input" tab
  showInputs: boolean = false;
  singleDfParam: boolean = false;
  // mappping of param to the switchces of their viewers/grids
  paramGridSwitches: Map<string, HTMLElement> = new Map();
  // mappping of param to their viewers
  paramViewers: Map<string, DG.Viewer[]> = new Map();
  // mappping of param to their html elements
  paramSpans: Map<string, HTMLElement> = new Map();
  // mappping of tab names to their content
  resultTabs: Map<String, HTMLElement> = new Map();
  resultsTabControl: DG.TabControl | undefined;
  resultsDiv: HTMLElement = ui.panel([], 'grok-func-results');
  controlsRoot: HTMLDivElement = ui.box();
  resultsRoot: HTMLDivElement = ui.box();
  historyRoot: HTMLDivElement = ui.divV([], {style: {'justify-content': 'center'}});
  inputsDiv: HTMLDivElement = ui.panel([], 'grok-func-results, ui-box');

  appendResultDataFrame(param: DG.FuncCallParam, options?: { caption?: string, category?: string, height?: number }) {
    const paramDf = param.value as DG.DataFrame;
    let caption = options?.caption;
    const height = options?.height ?? 400;
    const viewers: DG.Viewer[] = param.aux['viewers'] ?? []; // storing the viewers of Df
    caption ??= param.aux['viewerTitle'] ?? // description of the viewer ?
      ((this.singleDfParam && viewers.length == 1) ? '' : param.property.caption) ?? '';
    let existingViewers = this.paramViewers.get(param.name);
    if (existingViewers != null) {
      for (const v of existingViewers)
        v.dataFrame = paramDf;
      return;
    }
    existingViewers ??= [];
    this.paramViewers.set(param.name, existingViewers);
    if (viewers.length == 0)
      viewers.push(paramDf.plot.grid());

    const viewersToHide: HTMLElement[] = [];
    let blocks: number = 0;
    let blockWidth: number = 0;
    const gridWrapper = ui.box(null, {style: {height: '100%'}});

    const isGridSwitchable = wu(viewers).every((v: DG.Viewer) => v.type !== 'Grid');
    $(gridWrapper).hide();

    const getHeader = () => {
      const headerLabel = caption ?? paramDf.name ?? '';
      const header = ui.div([], 'grok-func-results-header');
      if (headerLabel != '') {
        const h = ui.h1(headerLabel);
        ui.Tooltip.bind(h, () => headerLabel);
        header.appendChild(h);
      }
      if (isGridSwitchable) {
        const icon = ui.iconSvg('table', (e) => {
          e.stopPropagation();
          ui.setDisplayAll(viewersToHide, viewersToHide.some((viewer) => viewer.style.display === 'none'));
          ui.setDisplay(gridWrapper, gridWrapper.style.display === 'none');
        }, 'Show grid');
        header.appendChild(icon);
        this.paramGridSwitches.set(param.name, icon);
      }
      if (!grok.shell.tables.includes(paramDf)) {
        header.appendChild(
          ui.icons.add((e: any) => {
            e.stopPropagation();
            const v = grok.shell.addTableView(paramDf);
            (async () => {
              for (const viewer of viewers) {
                if (viewer.type != 'Grid') {
                  const newViewer = await paramDf.plot.fromType(viewer.type) as DG.Viewer;
                  newViewer.setOptions(viewer.getOptions());
                  v.addViewer(newViewer);
                }
              }
            })();
          },
          'Add to workspace'
          ));
      }
      return header;
    };

    if (isGridSwitchable) {
      const grid = paramDf.plot.grid();
      gridWrapper.appendChild(grid.root);
      existingViewers.push(grid);
    }

    const header = getHeader();
    const wrapper = ui.divH([], {style: {height: '100%'}});
    viewers.forEach((viewer) => {
      if (viewer?.tags == null)
        return;
      existingViewers!.push(viewer);
      const viewerSize = viewer.tags['.block-size'] ?? 100; // getting `(block: 25)` value; if not defined, it is 100

      viewer.root.style.height = '100%';
      viewer.root.classList.add(`ui-block-${viewerSize}`);
      if (blocks + viewerSize <= 100) {
        viewersToHide.push(wrapper);
        blockWidth += viewerSize;
      }
      blocks += viewerSize;
      wrapper.appendChild(viewer.root);
    });

    const block = ui.divV([header, wrapper, gridWrapper], {style: {height: `${height}px`, width: `${blockWidth}%`}});
    block.classList.add('ui-box');
    this._appendResultElement(block, options?.category);
  }

  appendResultScalar(param: DG.FuncCallParam, options: { caption?: string, category: string, height?: number }) {
    const span = ui.span([`${options.caption ?? param.name}: `, `${param.value}`]);
    if (!this.paramSpans.get(param.name)) {
      if (this.resultsTabControl)
        this.resultTabs.get(options.category)!.appendChild(span);
      else
        this.resultsDiv.appendChild(span);
    } else {
      this.paramSpans.get(param.name)?.replaceWith(span);
    }
    this.paramSpans.set(param.name, span);
  }

  _appendResultElement(d: HTMLElement, category?: string) {
    if (category != null && this.resultsTabControl != undefined && this.resultTabs.get(category) != null)
      this.resultTabs.get(category)!.appendChild(d);
    else
      this.resultsDiv.appendChild(d);
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
