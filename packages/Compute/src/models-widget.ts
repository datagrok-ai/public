import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import $ from 'cash-dom';

export class ModelsWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    super(ui.panel());
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Model catalog');
    this.order = super.addProperty('order', DG.TYPE.STRING, '0');
    this.getModels();
  }

  async getModels(){
    let count = (await grok.dapi.scripts.filter('#model').list()).length; 
    let titles = [ui.h3('Simulation'), ui.h3('Other')]; 
    for (let i in titles){
        titles[i].style.color = 'var(--grey-3)';
        titles[i].style.fontSize = '12px';
        titles[i].style.textTransform = 'uppercase';
        titles[i].style.letterSpacing = '1px';
        titles[i].style.margin = '5px 0';
    }

    let modelsList = ui.divV([
        titles[0],
        ui.divV([],'models-widget-group'),
        titles[1],
        ui.divV([],'models-widget-group'),
    ]);

    grok.dapi.scripts
        .filter('#model')
        .list()
        .then((scripts) => scripts.map((p) => {
            if(p.tags.includes("simulation")){
                let div = ui.link('', p.path, undefined, null);
                div.append(modelData(p))
                $(modelsList).find('.models-widget-group')[0]?.append(div); 
            }else{               
                let div = ui.link('', p.path, undefined, null);
                div.append(modelData(p))
                $(modelsList).find('.models-widget-group')[1]?.append(div); 
            }
        }));

    this.root.appendChild(modelsList);
  }
}

function modelData(p: any){
    let img = ui.iconImage('logo', `/images/entities/${p.language}.png`);
    img.style.width = '30px';
    img.style.height = '30px';
    img.style.border = 'none';
    let card = ui.cards.summary(
        img,
        [
            ui.div(ui.link(p.friendlyName, p, undefined, null)),
            ui.div([p.updatedOn], {style:{color:'var(--grey-4)'}})
    ]);
    ui.tooltip.bind(card, p);
    return card;
  }