import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import wu from 'wu';
import '../css/model-card.css';


export class ModelCatalogView extends DG.CustomCardView {
  static findOrCreateCatalogView(
    viewName: string,
  ): ModelCatalogView {
    const modelsView = this.findModelCatalogView(viewName) ??
      this.createModelCatalogView(viewName);

    return modelsView as ModelCatalogView;
  }

  static createModelCatalogView(
    viewName: string,
  ): ModelCatalogView {
    const view = new this(viewName);
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
    viewName: string,
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

  async initRibbon() {
    // not used rn
  }

  async initMenu() {
    const standardHelpItems: {name: string, link: string}[] = [
      {
        name: 'Compute Engine',
        link: 'https://github.com/datagrok-ai/public/tree/master/packages/Compute',
      }, {
        name: 'Developing Models',
        link: 'https://datagrok.ai/help/compute/scripting'
      }
    ];
    const  customHelpItems = await this.getClientSpecificHelp();
    let help = this.ribbonMenu.group('Help');
    for (const item of [...standardHelpItems, ...customHelpItems])
      help = help.item(item.name,  () => window.open(item.link, '_blank'));
    help.endGroup();
  }

  private async getClientSpecificHelp(): Promise<{name: string, link: string}[]> {
    // use private ModelHub package
    const f: DG.Func = DG.Func.find({package: 'ModelHub', name: 'getHelpItems'})?.[0];
    if (!f)
      return [];
    const fc = f.prepare();
    await fc.call();
    return fc.getOutputParamValue() ?? [];
  }
}
