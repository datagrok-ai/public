import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {onboardModel} from './onboard-model';

/* eslint-disable */

export class ModelCatalogView extends DG.CustomCardView {

  constructor() {
    super({dataSource: grok.dapi.functions, permanentFilter: '#model'});

    this.root.classList.add('model-catalog-view');
    this.meta = new ModelHandler();
    this.name = 'Models';
    this.permanentFilter = '#model';
    this.renderMode = DG.RENDER_MODE.BRIEF;

    this.objectType = 'Func';
    this.categoryFilters = {
      'options.department': 'Department',
      'options.group': 'Group',
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
      'options.group',
      'options.process'
    ];
    this.showTree = true;
    this.initRibbon();
    this.initMenu();
    this.currentCall = grok.functions.getCurrentCall();
    //this.initFavorites().then((_) => {});
  }

  currentCall: DG.FuncCall;

  async initFavorites() {
    let fav = await grok.dapi.functions.filter('starredBy = @current #model').list();
    for (let f of fav) {
      //console.log(f.name);
      let fc = f.prepare();
      fc.parentCall = this.currentCall;
      let v = DG.View.create();
      v.name = f.name;
      v.parentView = this;
      v.parentCall = fc;

      grok.shell.addView(v);
      // @ts-ignore
      v.temp.model = f;
      ModelHandler.openModel(f, this.currentCall);
    }
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

  setCurrentModel(model: DG.Script) {
    let modelView = DG.View.forObject(model)!;
    let host: HTMLDivElement = this.root.querySelector('.grok-gallery-grid')! as HTMLDivElement;
    ui.empty(host);
    host.appendChild(modelView.root);
  }
}
