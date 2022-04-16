import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {onboardModel} from './onboard-model';

/* eslint-disable */

export class ModelCatalogView extends DG.CustomCardView {

  constructor() {
    super({dataSource: grok.dapi.functions, permanentFilter: '#model'});

    this.meta = new ModelHandler();
    this.name = 'Models';
    this.permanentFilter = '#model';

    this.objectType = 'Script';
    this.categoryFilters =       {
      'language': 'Language',
      'tag': 'Tags',
      'createdOn': 'Created',
      'author.id': 'Author'
    };
    this.initRibbon();
    this.initMenu();
  }

  initRibbon() {
    this.setRibbonPanels([[ui.icons.add(() => onboardModel())]]);
  }

  initMenu() {
    this.ribbonMenu
      .group('Models')
        .item('Add New...', () => onboardModel())
      .endGroup()
      .group('Help')
        .item('Compute Engine', () => window.open('https://github.com/datagrok-ai/public/tree/master/packages/Compute', '_blank'))
        .item('Developing Models', () => window.open('https://datagrok.ai/help/compute/scripting', '_blank'))
      .endGroup();
  }

  setCurrentModel(model: DG.Script) {
    let modelView = DG.View.forObject(model)!;
    let host: HTMLDivElement = this.root.querySelector('.grok-gallery-grid')! as HTMLDivElement;
    ui.empty(host);
    host.appendChild(modelView.root);
  }
}
