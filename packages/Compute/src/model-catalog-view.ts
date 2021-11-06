import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {ModelHandler} from "./model_handler";
let api = <any>window;

export class ModelCatalogView extends DG.CardView {

  constructor() {
    super(api.grok_CardView_Create({dataSource: grok.dapi.scripts, permanentFilter: '#model'}));

    this.meta = new ModelHandler();
    this.name = 'Models';
    this.permanentFilter = '#model';

    this.initToolbox();
    this.initRibbon();
  }

  initRibbon() {
    this.ribbonMenu
      .group('Models')
      .item('Add New...', () => grok.shell.info('Adding new model'));
  }

  initToolbox() {
    grok.dapi.scripts
      .filter('#model')
      .list()
      .then((models) => {
        this.toolbox = ui.divV(models.map(model => ui.render(model, {onClick: (_) => this.setCurrentModel(model)})));
      });
  }

  setCurrentModel(model: DG.Script) {
    let modelView = DG.View.forObject(model)!;
    let host: HTMLDivElement = this.root.querySelector('.grok-gallery-grid')!
    ui.empty(host);
    host.appendChild(modelView.root);
  }
}