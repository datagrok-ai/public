import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from "rxjs";
import {ModelHandler} from './model-handler';

/* eslint-disable */

let propsSub: rxjs.Subscription;

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
    this.currentCall = grok.functions.getCurrentCall();
    if (propsSub == null && window.localStorage.getItem('ModelCatalogShowProperties') == null) { // @ts-ignore
        propsSub = grok.events.onCurrentObjectChanged.subscribe((o) => {
                if (new ModelHandler().isApplicable(o.sender)) {
                  grok.shell.windows.showProperties = true;
                  window.localStorage.setItem('ModelCatalogShowProperties', 'false');
                  propsSub.unsubscribe();
                }
              })
      }
  }

  currentCall: DG.FuncCall;

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
