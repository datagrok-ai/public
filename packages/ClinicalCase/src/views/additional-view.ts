import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { getUniqueValues } from '../data-preparation/utils';


export class AdditionalView extends DG.ViewBase {

  tablesByArm: any;

  constructor(name) {
    super(name);
    this.name = name;
    
    let viewerTitle = {style:{
      'color':'var(--grey-6)',
      'margin':'12px 0px 6px 12px',
      'font-size':'16px',
    }};

    const plot = DG.Viewer.boxPlot(study.domains.dm, {
      category: 'ARM',
      value: 'AGE',
    });
    plot.root.prepend(ui.divText('Age distribution by treatment arm', viewerTitle));
   
    this.root.appendChild(plot.root);

  } 
}