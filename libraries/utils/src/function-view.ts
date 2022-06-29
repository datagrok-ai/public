/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import wu from 'wu';
import $ from 'cash-dom';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class FunctionView extends DG.ViewBase {
  constructor(funcCall?: DG.FuncCall) {
    super();
    this.box = true;
    this.context = DG.Context.cloneDefault();

    if (!funcCall) return;

    this.linkFunccall(funcCall);
    this.init();
    this.build();
  }

  get func() {
    return this.funcCall?.func;
  }

  private _type: string = 'function';
  public get type(): string {
    return this._type;
  }

  public get isInputPanelRequired() {
    return this.func?.inputs.some((p) => p.propertyType == DG.TYPE.DATA_FRAME && p.options['viewer'] != null) || false;
  }

  public get outParamCategories() {
    return [
      ...new Set(this.func!.outputs.map((p) => p.category)) // get all output params' categories
    ]; // keep only unique of them
  }

  public get outputTabsLabels() {
    return [
      ...this.outParamCategories,
      ...this.outParamCategories.find((val) => val === 'Misc') ? ['Output'] : [], // if no categories are stated, the default category is added
    ];
  }

  public get tabsLabels() {
    return Object.keys(this.paramToCategoryMap);
  }

  public get paramToCategoryMap() {
    const map = {} as Record<string, string[]>;
    if (this.isInputPanelRequired)
      this.func!.inputs.forEach((p) => map['Input'] ? map['Input'].push(p.name): map['Input'] = [p.name]);

    this.func!.outputs.forEach((p) => map[p.category === 'Misc' ? 'Output': p.category] ? map[p.category === 'Misc' ? 'Output': p.category].push(p.name) : map[p.category === 'Misc' ? 'Output': p.category] = [p.name]);

    return map;
  }

  readonly context: DG.Context;
  protected _funcCall?: DG.FuncCall;
  protected lastCall?: DG.FuncCall;

  /** Link FuncCall to the view */
  public linkFunccall(funcCall: DG.FuncCall) {
    this._funcCall = funcCall;
  }

  /** Get copy of FuncCall of the view */
  public get funcCall(): DG.FuncCall | undefined {
    return this._funcCall;
  }

  /** Init custom logic */
  init() {
    if (this.funcCall && this.func) {
      this.funcCall.aux['view'] = this;
      this.funcCall.context = this.context;

      for (const inParam of wu(this.funcCall.inputParams.values() as DG.FuncCallParam[])
        .filter((p: DG.FuncCallParam) =>
          p.property.propertyType == DG.TYPE.DATA_FRAME && p.property.options['viewer'] != null)) {
        this.subs.push(inParam.onChanged.subscribe(async (param: DG.FuncCallParam) => {
          this.clearResults(false);
          param.processOutput();

          this.appendOutputDf(param, {
            height: 400,
            category: 'Input'
          });
        }));
      }
    }
  }

  /** Override to create a fully custom UI */
  build(): void {
    this.root.appendChild(this.buildIO());

    this.buildRibbonPanels();
    this.buildHistoryBlock();
    this.buildRibbonMenu();
  }

  /** Override to create a custom input-output block */
  buildIO(): HTMLElement {
    const inputBlock = this.buildInputBlock();
    const outputBlock = this.buildOutputBlock();

    return ui.splitH([inputBlock, outputBlock]);
  }

  /** Override to create a custom input block. */
  buildInputBlock(): HTMLElement {
    if (!this.funcCall) return this.controlsRoot;

    const funcDiv = ui.div([this.renderRunSection(this.funcCall)], 'ui-div');
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
    return this.getRibbonPanels();
  }

  /** Override to create a custom ribbon menu on the top. */
  buildRibbonMenu() {

  }

  /** Saves the computation results to the historical results, returns its id. See also {@link loadRun}. */
  async saveRun(): Promise<DG.FuncCall> {
    //@ts-ignore
    return await grok.dapi.functions.calls.save(this.funcCall!);
  }

  /** Loads the specified historical results. See also {@link saveRun}. */
  async loadRun(funcCallId: string): Promise<DG.FuncCall> {
    return await grok.dapi.functions.calls.include('inputs, outputs').find(funcCallId);
  }

  /** Loads all the function call of this function. */
  /** ACHTUNG: FuncCall inputs/outputs are not included */
  async pullRuns(): Promise<DG.FuncCall[]> {
    return await grok.dapi.functions.calls.filter(`func.id="${this.func?.id}"`).list();
  }

  clearResults(switchToOutput: boolean = true) {
    const categories = this.tabsLabels;

    if ((categories.length > 1 || (categories.length == 1 && categories[0] != 'Misc'))) {
      this.resultsTabControl = DG.TabControl.create();
      for (const c of categories) {
        if (!this.resultTabs.has(c))
          this.resultTabs.set(c, ui.div([], 'ui-panel, grok-func-results, ui-box'));
        let name = c;
        if (this.isInputPanelRequired && categories.length == 2 && c == 'Misc')
          name = 'OUTPUT';
        this.resultsTabControl.addPane(name, () => this.resultTabs.get(c) ?? ui.div());
      }
      if (categories.length > 1 && this.isInputPanelRequired && switchToOutput)
        this.resultsTabControl.currentPane = this.resultsTabControl.panes[1];
      this.resultsDiv = this.resultsTabControl.root;
    }
  }

  async run(): Promise<void> {
    if (!this.funcCall) throw new Error('The correspoding function is not specified');

    try {
      await this.compute(this.funcCall);
    } catch (e) {

    }
    this.lastCall = this.funcCall;

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
        this.appendOutputDf(p, {caption: p.property.name, category: p.property.category});
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

  // mappping of param to the switchces of their viewers/grids
  paramGridSwitches: Map<string, HTMLElement> = new Map();
  // mappping of param to their viewers placed on DOM
  existingParamViewers: Map<string, DG.Viewer[]> = new Map();
  // mappping of param to their html elements
  paramSpans: Map<string, HTMLElement> = new Map();
  // mappping of tab names to their content
  resultTabs: Map<String, HTMLElement> = new Map();
  resultsTabControl: DG.TabControl | undefined;
  resultsDiv: HTMLElement = ui.panel([], 'grok-func-results');
  controlsRoot: HTMLDivElement = ui.box(null, {style: {maxWidth: '370px'}});
  resultsRoot: HTMLDivElement = ui.box();
  historyRoot: HTMLDivElement = ui.divV([], {style: {'justify-content': 'center'}});
  inputsRoot: HTMLDivElement = ui.panel([], 'grok-func-results, ui-box');

  appendOutputDf(param: DG.FuncCallParam, options?: { caption?: string, category?: string, height?: number }) {
    const paramDf = param.value as DG.DataFrame;
    const height = options?.height ?? 400;
    const paramViewers: DG.Viewer[] = param.aux['viewers'] ?? []; // storing the viewers of Df
    const paramCaption = options?.caption ?? param.aux['viewerTitle'] ?? param.property.caption ?? '';
    let existingViewers = this.existingParamViewers.get(param.name);
    if (existingViewers != null) {
      for (const v of existingViewers)
        v.dataFrame = paramDf;
      return;
    }
    existingViewers ??= [];
    this.existingParamViewers.set(param.name, existingViewers);
    if (paramViewers.length == 0)
      paramViewers.push(paramDf.plot.grid());

    const viewersToHide: HTMLElement[] = [];
    let blocks: number = 0;
    let blockWidth: number = 0;
    const gridWrapper = ui.box(null, {style: {height: '100%'}});

    const isGridSwitchable = wu(paramViewers).every((v: DG.Viewer) => v.type !== 'Grid');
    $(gridWrapper).hide();

    const getHeader = () => {
      const headerLabel = paramCaption ?? paramDf.name ?? '';
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
              for (const viewer of paramViewers) {
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
    paramViewers.forEach((viewer) => {
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
