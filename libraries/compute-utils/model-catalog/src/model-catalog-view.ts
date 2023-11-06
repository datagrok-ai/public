import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {ModelHandler} from './model-handler';
import wu from 'wu';

/* eslint-disable */

let propsSub: rxjs.Subscription;

export class ModelCatalogView extends DG.CustomCardView {
  static findOrCreateCatalogView(
    viewName: string,
    funcName: string,
    currentPackage: DG.Package,
  ): ModelCatalogView {
    let modelsView = this.findModelCatalogView(viewName);
    if (modelsView == null) {
      modelsView = this.createModelCatalogView(viewName, funcName, currentPackage)
    }
    return modelsView as ModelCatalogView;
  }

  private static createModelCatalogView(viewName: string, funcName: string, currentPackage: DG.Package): ModelCatalogView {
    const mc: DG.Func = DG.Func.find({package: currentPackage.name, name: viewName})[0];
    const fc = mc.prepare();
    const view = new this(viewName, funcName, currentPackage);
    view.parentCall = fc;
    return view;
  }

  private static findModelCatalogView(
    viewName: string,
  ): ModelCatalogView | undefined {
    return wu(grok.shell.views).find((v) => v.parentCall?.func.name == viewName) as ModelCatalogView;
  }

  static openWithModel(
    viewName: string,
    funcName: string,
    currentPackage: DG.Package,
  ) {
    const modelsView = this.findOrCreateCatalogView(viewName, funcName, currentPackage);
  
    const pathSegments = window.location.pathname.split('/');
    if (pathSegments.length > 3) {
      grok.dapi.functions.filter(`shortName = "${pathSegments[3]}" and #model`).list().then((lst) => {
        if (lst.length == 1)
          ModelHandler.openModel(lst[0]);
      });
    } else {
      grok.shell.v = modelsView;
    }
  }

  private bindModel(fc: DG.FuncCall) {
    const modelsView = ModelCatalogView.findOrCreateCatalogView(this.viewName, this.funcName, this.currentPackage);
    fc.aux['showOnTaskBar'] = false;
    if (modelsView != null) {
      const parentCall = modelsView.parentCall;
      if (parentCall != null) {
        parentCall.aux['view'] = modelsView;
        fc.parentCall = parentCall;
      }
    }
  }

  constructor(
    private viewName: string,
    private funcName: string,
    private currentPackage: DG.Package,
  ) {
    super({dataSource: grok.dapi.functions, permanentFilter: '#model'});

    this.root.classList.add('model-catalog-view');
    this.name = viewName;
    this.permanentFilter = '#model';
    this.renderMode = DG.RENDER_MODE.BRIEF;

    this.objectType = 'Func';
    this.categoryFilters = {
      'options.department': 'Department',
      'options.HL_process': 'HL_process',
      'options.process': 'Process',
      'tag': 'Tags'
    };

    this.filters = {
      'All': '',
      'Favorites': 'starredBy = @current',
      'Used by me': 'usedBy = @current'
    };

    this.hierarchy = [
      'options.department',
      'options.HL_process',
      'options.process'
    ];

    this.hierarchyProperties = {
      'options.department': 'Department',
      'options.HL_process': 'HL_process',
      'options.process': 'Process'
    };

    this.showTree = true;
    this.initRibbon();
    this.initMenu();
    if (propsSub == null && window.localStorage.getItem('ModelCatalogShowProperties') == null) { // @ts-ignore
        propsSub = grok.events.onCurrentObjectChanged.subscribe((o) => {
          if (ModelHandler.isApplicable(o.sender)) {
            grok.shell.windows.showProperties = true;
            window.localStorage.setItem('ModelCatalogShowProperties', 'false');
            propsSub.unsubscribe();
          }
        })
      }
    grok.shell.windows.showHelp = false;

    setTimeout(async () => {
      grok.functions.onBeforeRunAction.subscribe((fc) => {
        if (fc.func.hasTag('model') || fc.func.hasTag('model-editor'))
          this.bindModel(fc);
      });
    });
  }

  initRibbon() {
    this.setRibbonPanels([[ui.icons.sync(() => this.refresh())]]);
  }

  initMenu() {
    this.ribbonMenu
      .group('Help')
        .item('Compute Engine', () => window.open('https://github.com/datagrok-ai/public/tree/master/packages/Compute', '_blank'))
        .item('Developing Models', () => window.open('https://datagrok.ai/help/compute/scripting', '_blank'))
      .endGroup();
  }
}
