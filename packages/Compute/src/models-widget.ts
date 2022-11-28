import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {ModelHandler} from './model-handler';
import $ from 'cash-dom';
import {FunctionView} from "@datagrok-libraries/utils/src/function-view";
const api = <any>window;

/*

export class ModelsWidget extends DG.Widget {
  caption: string;
  order: string;
  meta: ModelHandler;
  name: string;
  permanentFilter: string;
  objectType: string;
  modelsList: HTMLDivElement;

  constructor() {
    super(ui.box());
    
    this.meta = new ModelHandler();
    this.name = 'Models';
    this.permanentFilter = '#model';
    this.objectType = 'Script';

    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Models');
    this.order = super.addProperty('order', DG.TYPE.STRING, '0');

    let topbar = ui.divH([
        ui.choiceInput('', 'Recent', ['Recent','Newest','Oldest'], (v:string) => {
            this.modelsList.innerHTML = '';
            this.getModels(v);
        }).root,
        ui.link('Open Model Catalog', ()=>{
            grok.functions.call("Compute:ModelCatalog");
            setTimeout(() => {
            for (let v of grok.shell.views)
                if (v.name === 'Models')
                    grok.shell.v = v;
            }, 500);
    })
    ]);
    topbar.style.justifyContent = 'space-between';
    topbar.style.margin = '0 12px';
    topbar.style.alignItems = 'center';

    this.modelsList = ui.divV([]);
    this.modelsList.style.marginLeft = '12px';

    this.root.appendChild(ui.splitV([
        ui.box(topbar, {style:{maxHeight:'40px'}}),
        this.modelsList
    ]));

    DG.ObjectHandler.register(new ModelHandler());

    this.getModels('updatedOn');
  }

    

  async getModels(order:string){
    let b = true;
    let v = order;

    switch(order){
        case 'Recent': 
            v = 'updatedOn';
            break;
        case 'Newest': 
            v = 'createdOn';
            break;
        case 'Oldest': 
            v = 'createdOn';
            b = false;
            break;
        default:
            v = 'updatedOn';     
    }

    grok.dapi.scripts
        .filter('#model')
        .order(v, b)
        .list()
        .then((scripts) => scripts.map((p) => {
            let d = this.meta.renderMarkup(p);
            d.ondblclick = (e) => {
                ModelsWidget.openModel(p);
              }
            d.style.color = 'var(--blue-1)';
            d.style.margin = '4px 0';
            $(d).find('.ui-label').addClass('ui-link');
            
            this.modelsList.append(d)  
        }));
  }

  static openModel(x: DG.Script, parentCall?: DG.FuncCall) {
    let views = []  
    for (let v of grok.shell.views)
      views.push(v.name)
    
    if (views.includes('Models')){
          let view = new FunctionView(x);
          view.parentCall = parentCall!;
          grok.shell.addView(view);
    }else{
          grok.functions.call("Compute:ModelCatalog");
          setTimeout(() => {
          let view = new FunctionView(x);
          view.parentCall = parentCall!;
          grok.shell.addView(view);
          }, 500)
    }
  }

} */
