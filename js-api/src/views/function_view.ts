import {Func, Property} from "../entities";
import {toJs} from "../wrappers";
import {Context, FuncCall, FuncCallParam} from "../functions";
import * as ui from "../../ui";
import {Viewer} from "../viewer";
import {TabControl} from "../widgets";
import {div, iconSvg, setDisplay, setDisplayAll, splitH, Tooltip} from "../../ui";
import {DataFrame} from "../dataframe";
import {View, ViewBase} from "./view";
import {TYPE} from "../const";
import wu from "wu";
import $ from "cash-dom";
import {shell} from "../../grok";


let api = <any>window;

export class FunctionView extends ViewBase {
  constructor(func: Func) {
    super();
    this.func = func;
    this.context = Context.create();
    this.controlsRoot.style.maxWidth = '370px';
    this.box = true;
    this._setFunction(undefined).then(r => {});
  }
  public func: Func;
  private readonly context: Context;


  async _setFunction(call: FuncCall | undefined) {

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

    this.singleDfParam = wu(this.func.outputs).filter((p) => p.propertyType == TYPE.DATA_FRAME).toArray().length == 1;
    let runSection = this.renderRunSection(call);
    //var fullRunSection = div();

    let funcDiv = div([runSection], 'ui-div');
    //var fullDiv = div('ui-panel', [getFullSection(func, call, fullRunSection)]);
    //htmlSetDisplay(fullDiv, false);

     for (let inParam of wu(call.inputParams.values() as FuncCallParam[]).filter((p: FuncCallParam) => p.property.propertyType == TYPE.DATA_FRAME && p.property.options['viewer'] != null)) {
       this.showInputs = true;
       let self = this;
       this.subs.push(inParam.onChanged.subscribe(async function(param: FuncCallParam) {
         self.clearResults(true, false);
         param.processOutput();

         self.appendResultDataFrame(param,
           { height: (self.singleDfParam && !shell.windows.presentationMode) ? 600 : 400,
           category: 'Input'});
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

    this.root.appendChild(splitH([this.controlsRoot, this.resultsRoot], ));
   /*
    Routing.setPath(this.path, f.name, true);
    */
    this.resultsRoot.innerHTML = '';
    this.resultsRoot.appendChild(this.resultsDiv);
    this.name = this.func.friendlyName;
  }

  clearResults(clearTabs: boolean = true, showOutput: boolean = true) {
    let categories: string[] = [];
    //  resultTabs.clear();
    if (this.showInputs) {
      categories.push('Input');
      this.resultTabs.set('Input', this.inputsDiv);
    }
    for (let p of this.func.outputs) {
      if (categories.includes(p.category))
        continue;
      categories.push(p.category);
    }
    if (clearTabs && (categories.length > 1 || (categories.length == 1 && categories[0] != 'Misc'))) {
      this.resultsTabControl = TabControl.create();
      for (let c of categories) {
        if (!this.resultTabs.has(c))
          this.resultTabs.set(c, div([],'ui-panel,grok-func-results'));
        let name = c;
        if (this.showInputs && categories.length == 2 && c == 'Misc')
          name = 'Output';
        this.resultsTabControl.addPane(name, () => this.resultTabs.get(c) ?? div());
      }
      if (categories.length > 1 && this.showInputs && showOutput)
        this.resultsTabControl.currentPane = this.resultsTabControl.panes[1];
      this.resultsDiv = this.resultsTabControl.root;
    }
    else {
      this.resultsDiv = div([], 'ui-panel,grok-func-results');
      this.resultsTabControl = undefined;
    }
    this.resultsRoot.innerHTML = '';
    this.resultsRoot.appendChild(this.resultsDiv);
  }

  renderRunSection(call: FuncCall): HTMLElement {
    return ui.wait(async () => {
      let runButton = ui.bigButton('Run', async () => {
        call.aux['view'] = this.dart;
        await call.call(true, undefined, {processed: true});
        console.log(call);
        for (let [name, p] of call.outputParams) {
          console.log(p);
          p.processOutput();

          if (p.property.propertyType == TYPE.DATA_FRAME) {
            console.log('add df');
          this.appendResultDataFrame(p, { caption: p.property.name, category: p.property.category});
          }
        }
      });
      let editor: HTMLDivElement = await call.getEditor(true);
      editor.appendChild(ui.buttonsInput([runButton]));
      return editor;
    });
  }

  showInputs: boolean = false;
  singleDfParam: boolean = false;
  paramViewers: Map<string, Viewer[]> = new Map();
  resultsTabControl: TabControl | undefined;
  resultTabs: Map<String, HTMLElement> = new Map();
  resultsDiv: HTMLElement = div([],'ui-panel,grok-func-results');
  controlsRoot: HTMLDivElement = div([], 'ui-box');
  resultsRoot: HTMLDivElement = div([], 'ui-box');
  inputsDiv: HTMLDivElement = div([], 'ui-panel,grok-func-results');

  appendResultDataFrame(param: FuncCallParam, options?: { caption?: string, category?: string, height?: number}) {
    let df = param.value;
    let caption = options?.caption;
    let height = options?.height ?? 400;
    let viewers: Viewer[] = param.aux['viewers'] ?? [];
    caption ??= param.aux['viewerTitle'] ?? ((this.singleDfParam && viewers.length == 1) ? '' : param.property.caption) ?? '';
    let existingViewers: Viewer[] | undefined = this.paramViewers.get(param.name);
    if (existingViewers != null) {
      for (let v of existingViewers) {
        v.dataFrame = df;
      }
      return;
    }
    existingViewers ??= [];
    this.paramViewers.set(param.name, existingViewers);
    if (viewers.length == 0)
      viewers.push(df.plot.grid());

    let hideList: HTMLElement[] = [];
    let blocks: number = 0;
    let blockSize: number = 0;
    let gridWrapper = div([], 'ui-div,ui-block');

    console.log(viewers);
    let gridSwitch = !wu(viewers).some((v: Viewer) => v.type == 'grid');
    $(gridWrapper).hide();

    let getHeader = (sw: boolean) => {
      let s = caption ?? df.name ?? '';
      let header = div([], 'grok-func-results-header');
      if (s != '') {
        let h = ui.h1(s);
        Tooltip.bind(h, () => s);
        header.appendChild(h);
      }
      if (gridSwitch) {
        let icon = iconSvg('table', (e) => {
            e.stopPropagation();
            setDisplayAll(hideList, !sw);
            setDisplay(gridWrapper, sw);
            gridWrapper.classList.add('ui-block-$blockSize');
          }
          , 'Show grid');
        if (!sw)
          icon.classList.add('active');
        header.appendChild(icon);
      }
      if (!shell.tables.includes(df))
        header.appendChild(ui.icons.add((e: any) => {
            e.stopPropagation();
            let v = shell.addTableView(df);
            (async () => {
              for (let viewer of viewers) {
                if (viewer.type != 'grid') {
                  let newViewer = await df.plot.fromType(viewer.type) as Viewer;
                  newViewer.setOptions(viewer.getOptions());
                  v.addViewer(newViewer);
                }
              }
            })();
          },
          'Add to workspace'
        ));
      return header;
    }

    if (gridSwitch) {
      gridWrapper.appendChild(getHeader(false));
      let grid = df.plot.grid();
      grid.root.style.height = `${height}px`;
      gridWrapper.appendChild(grid.root);
      existingViewers.push(grid);
    }

    let header = getHeader(true);
    this._appendResultElement(gridWrapper, options?.category);
    for (let viewer of viewers) {
      existingViewers.push(viewer);
      let block = viewer.tags['.block-size'] ?? 100;
      let wrapper = div([], `ui-div,ui-block,ui-block-${block}`);
      if (blocks + block <= 100) {
        hideList.push(wrapper);
        blockSize += block;
      }
      blocks += block;
      wrapper.appendChild(header);
      header = div([], 'grok-func-results-header');
      wrapper.appendChild(viewer.root);
      if (viewer.type == 'grid') {
        // @ts-ignore
        let totalHeight = viewer.getOptions()['rowHeight'] *
          (viewer.dataFrame!.rowCount + 1);
        if (totalHeight < height)
          height = totalHeight;
      }
      viewer.root.style.height = `${height}px`;
      this._appendResultElement(wrapper, options?.category);
      console.log('add viewer');
    }

  }

  _appendResultElement(d: HTMLElement, category?: string) {
    if (category != null && this.resultsTabControl != undefined && this.resultTabs.get(category) != null)
      this.resultTabs.get(category)!.appendChild(d);
    else
      this.resultsDiv.appendChild(d);
  }
}
