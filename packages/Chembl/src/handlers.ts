import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SemanticValue, x} from "datagrok-api/dg";

export class ChemblIdHandler extends DG.ObjectHandler {
  get type(): string { return 'CHEMBL_ID'; }

  isApplicable(x: any): boolean {
    return x instanceof DG.SemanticValue && x.semType == 'CHEMBL_ID';
  }

  renderInnerCard(x: any) {
    const id = (x as SemanticValue).value;
    return ui.divV([
      ui.h3(id),
      ui.wait(async () => ui.bind(x, grok.chem.drawMolecule(await grok.functions.call('Chembl:chemblIdToSmiles', {id: id}))))
    ], { style: { width: '220px', height: '150px' }});
  }

  renderCard(x: any) {
    return ui.card(this.renderInnerCard(x));
  }

  renderTooltip(x: any, context?: any): HTMLElement {
    return this.renderInnerCard(x);
  }

  renderProperties(x: any, context?: any): HTMLElement {
    return ui.panels.infoPanel(x).root;
  }
}