import wu from 'wu';
import $ from 'cash-dom';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class FunctionView extends DG.ViewBase {
  constructor(func: DG.Func) {
    super();
    this.func = func;
    this.context = DG.Context.cloneDefault();
    this.controlsRoot.style.maxWidth = '370px';
    this.box = true;
    this._setFunction(undefined).then((r) => {});
  }

  private _type: string = 'function';
  public get type(): string {
    return this._type;
  }
  public func: DG.Func;
  private readonly context: DG.Context;


  async _setFunction(call: DG.FuncCall | undefined) {
    /*
      var meta = EntityMeta.forEntity(f);
      func = await meta.refresh(f);
    */
    if (call != undefined) {
      /*  if (call.validateParameterValues().isEmpty)
          runFunc(call).then((_) => call.reset());
        _callPath = '/${call.toUri()}';
        Routing.setPath(this.path, f.name, true);*/
    }

    call ??= this.func.prepare();
    call.aux['view'] = this;
    call.context = this.context;

    this.singleDfParam = wu(this.func.outputs).filter((p) => p.propertyType == DG.TYPE.DATA_FRAME)
      .toArray().length == 1;
    const runSection = this.renderRunSection(call);
    //var fullRunSection = div();

    const funcDiv = ui.div([runSection], 'ui-div');
    //var fullDiv = div('ui-panel', [getFullSection(func, call, fullRunSection)]);
    //htmlSetDisplay(fullDiv, false);

    for (const inParam of wu(call.inputParams.values() as DG.FuncCallParam[])
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
    }

    this.controlsRoot.innerHTML = '';
    this.controlsRoot.appendChild(funcDiv);
    //this.controlsRoot.appendChild(fullDiv);

    /*advancedModeSwitch.onChanged.listen((_) {
      if (advancedModeSwitch.value)
        fullRunSection.append(runSection);
      else
        funcDiv.append(runSection);
      htmlSetDisplay(fullDiv, advancedModeSwitch.value);
      htmlSetDisplay(funcDiv, !advancedModeSwitch.value);
    });*/
    //setRibbon();

    this.root.appendChild(ui.splitH([this.controlsRoot, this.resultsRoot] ));
    /*
    Routing.setPath(this.path, f.name, true);
    */
    this.resultsRoot.innerHTML = '';
    this.resultsRoot.appendChild(this.resultsDiv);
    this.name = this.func.friendlyName;
  }

  clearResults(clearTabs: boolean = true, showOutput: boolean = true) {
    const categories: string[] = [];
    //  resultTabs.clear();
    if (this.showInputs) {
      categories.push('INPUT');
      this.resultTabs.set('INPUT', this.inputsDiv);
    }
    for (const p of this.func.outputs) {
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
    this.resultsRoot.innerHTML = '';
    this.resultsRoot.appendChild(this.resultsDiv);
  }

  renderRunSection(call: DG.FuncCall): HTMLElement {
    return ui.wait(async () => {
      const runButton = ui.bigButton('Run', async () => {
        call.aux['view'] = this.dart;
        await call.call(true, undefined, {processed: true});
        this.clearResults(false, true);
        for (const [, p] of call.outputParams) {
          console.log('param', p.value);
          p.processOutput();
          if (p.property.propertyType == DG.TYPE.DATA_FRAME && p.value != null)
            this.appendResultDataFrame(p, {caption: p.property.name, category: p.property.category});
        }
      });
      const editor: HTMLDivElement = await call.getEditor(true);
      editor.appendChild(ui.buttonsInput([runButton]));
      return editor;
    });
  }

  showInputs: boolean = false;
  singleDfParam: boolean = false;
  paramViewers: Map<string, DG.Viewer[]> = new Map();
  resultsTabControl: DG.TabControl | undefined;
  resultTabs: Map<String, HTMLElement> = new Map();
  resultsDiv: HTMLElement = ui.panel([], 'grok-func-results');
  controlsRoot: HTMLDivElement = ui.box();
  resultsRoot: HTMLDivElement = ui.box();
  inputsDiv: HTMLDivElement = ui.panel([], 'grok-func-results');

  appendResultDataFrame(param: DG.FuncCallParam, options?: { caption?: string, category?: string, height?: number}) {
    const df = param.value;
    console.log('dataframe', df);
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

    const gridSwitch = !wu(viewers).some((v: DG.Viewer) => v.type == 'grid');
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
              if (viewer.type != 'grid') {
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
}
