/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TemplateHandler} from "../../HitTriage/src/template_handler";
import {ModelHandler} from "./model_handler";

export let _package = new DG.Package();

//tags: init
export function init() {
  console.log('init');
  DG.ObjectHandler.register(new ModelHandler());
}

//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}


//name: Model Catalog
//tags: app
export function modelCatalog() {
  let v = DG.CardView.create({dataSource: grok.dapi.scripts, permanentFilter: '#modelhub'})
  v.meta = new ModelHandler();
  v.name = 'Models';
  v.permanentFilter = '#modelhub';
  grok.shell.addView(v);
}