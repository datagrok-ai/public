import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import wu from 'wu';
import '../css/model-card.css';


export class ModelCatalogView extends DG.CustomCardView {
  static findOrCreateCatalogView(
    viewName: string,
    funcName: string,
    currentPackage: DG.Package,
  ): ModelCatalogView {
    const modelsView = this.findModelCatalogView(viewName) ??
      this.createModelCatalogView(viewName, funcName, currentPackage);

    return modelsView as ModelCatalogView;
  }

  static createModelCatalogView(
    viewName: string,
    funcName: string,
    currentPackage: DG.Package,
  ): ModelCatalogView {
    const mc: DG.Func = DG.Func.find({package: currentPackage.name, name: funcName})[0];
    const fc = mc.prepare();
    const view = new this(viewName);
    view.parentCall = fc;
    return view;
  }

  static findModelCatalogView(
    viewName: string,
  ): ModelCatalogView | undefined {
    return wu(grok.shell.views).find((v) => v.name == viewName) as ModelCatalogView;
  }

  static bindModel(view: DG.ViewBase, fc: DG.FuncCall) {
    fc.aux['showOnTaskBar'] = false;

    const parentCall = view.parentCall;
    if (parentCall != null) {
      parentCall.aux['view'] = view;
      fc.parentCall = parentCall;
    }
  }

  constructor(
    private viewName: string,
  ) {
    super({dataSource: grok.dapi.functions, permanentFilter: '#model'});

    this.root.classList.add('model-catalog-view');
    this.meta = new ModelHandler();
    this.name = viewName;
    this.permanentFilter = '#model';
    this.renderMode = DG.RENDER_MODE.BRIEF;

    // Smth is wrong with custom methods after getting an instance via grok.shell.views

    this.objectType = 'Func';
    this.categoryFilters = {
      'options.department': 'Department',
      'options.HL_process': 'HL_process',
      'options.process': 'Process',
      'tag': 'Tags',
    };

    this.filters = {
      'All': '',
      'Favorites': 'starredBy = @current',
      'Used by me': 'usedBy = @current',
    };

    this.hierarchy = [
      'options.department',
      'options.HL_process',
      'options.process',
    ];

    this.hierarchyProperties = {
      'options.department': 'Department',
      'options.HL_process': 'HL_process',
      'options.process': 'Process',
    };

    this.showTree = true;
    this.initRibbon();
    this.initMenu();
    grok.shell.windows.showProperties = false;
  }

  initRibbon() {
    // place additional icons here if necessary
  }

  initMenu() {
    this.ribbonMenu
      .group('Help')
      .item('Compute Engine',
        () => window.open('https://github.com/datagrok-ai/public/tree/master/packages/Compute', '_blank'))
      .item('Developing Models',
        () => window.open('https://datagrok.ai/help/compute/scripting', '_blank'))
      .endGroup();
  }
}
