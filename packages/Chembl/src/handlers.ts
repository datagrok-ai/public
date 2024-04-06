import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SemanticValue, x} from "datagrok-api/dg";

export class ChemblIdHandler extends DG.ObjectHandler {
  get type(): string { return 'CHEMBL_ID'; }

  isApplicable(x: any): boolean {
    return x instanceof DG.SemanticValue && x.semType == 'CHEMBL_ID';
  }

  renderCard(x: any) {
    const id = (x as SemanticValue).value;
    return ui.divV([
      id,
      ui.wait(async () => ui.bind(x, grok.chem.drawMolecule(await grok.functions.call('Chembl:chemblIdToSmiles', {id: id}))))
    ]);
  }

  renderTooltip(x: any, context?: any): HTMLElement {
    return this.renderCard(x);
  }

  renderProperties(x: any, context?: any): HTMLElement {
    return ui.panels.infoPanel(x).root;
  }
}