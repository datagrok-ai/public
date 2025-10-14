import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SemanticValue, x} from 'datagrok-api/dg';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {moleculeRenderer} from '@datagrok-libraries/db-explorer/src/renderer';

export class ChemblIdHandler extends DG.ObjectHandler {
  get type(): string {return 'CHEMBL_ID';}

  isApplicable(x: any): boolean {
    return x instanceof DG.SemanticValue && x.semType == 'CHEMBL_ID';
  }

  renderInnerCard(x: any) {
    const id = (x as SemanticValue).value;
    return ui.divV([
      ui.h3(id),
      ui.wait(async () =>
        ui.bind(x, grok.chem.drawMolecule(await grok.functions.call('Chembl:chemblIdToSmiles', {id: id})))),
    ], {style: {width: '220px', height: '150px'}});
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

export async function registerChemblIdHandler(_package: DG.Package) {
  // chech if chembl connection is available
  const chembl = await grok.dapi.connections.filter('name="CHEMBL"').first();
  if (!chembl) {
    console.warn('CHEMBL connection not found. CHEMBL object handlers not registered');
    return;
  }
  const exp = await DBExplorer.initFromConfigPath(_package);
  if (!exp) {
    grok.shell.error('Failed to load db-explorer config');
    return;
  }
  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'structure' || lc.includes('smiles')) && typeof value === 'string' && grok.chem.checkSmiles(value);
  }, (value) => moleculeRenderer(value as string));

  console.log('CHEMBL object handlers registered');
}
