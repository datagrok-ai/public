import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as crud from "../plates-crud";
import { propertySchemaView } from '../views/plates-schema-view';

export class PlateTemplateHandler extends DG.ObjectHandler<DG.SemanticValue<crud.PlateTemplate>> {
  
  get type() { return crud.TYPE.TEMPLATE; }
  get description() { return 'A set of properties for plates and wells'; }

  isApplicable(x: any): boolean {
    return x instanceof DG.SemanticValue && crud.TYPE.TEMPLATE === x.semType;
  }

  renderView(x: DG.SemanticValue<crud.PlateTemplate>, context: any): HTMLElement {
    return propertySchemaView(x.value).root;
  }

  renderTooltip(x: DG.SemanticValue<crud.PlateTemplate>, context: any): HTMLElement {
    return ui.div(x.value.description);
  }

  constructor() {
    super();

    this.registerParamFunc('Edit...',  (template: DG.SemanticValue<crud.PlateTemplate>) => PlateTemplateHandler.editTemplate(template.value));
    this.registerParamFunc('Clone...', (template: DG.SemanticValue<crud.PlateTemplate>) => PlateTemplateHandler.openCloneTemplateView(template.value));
    this.registerParamFunc('Delete', async (template: DG.SemanticValue<crud.PlateTemplate>) => await PlateTemplateHandler.tryDeleteTemplate(template.value));
  }

  static editTemplate(template: crud.PlateTemplate) {
    grok.shell.addPreview(propertySchemaView(template));
  }

  static async openCloneTemplateView(template: crud.PlateTemplate) {
    grok.shell.addPreview(propertySchemaView({...template, id: -1}));
  }

  static async tryDeleteTemplate(template: crud.PlateTemplate): Promise<boolean> {
    if (await crud.plateTemplatePropertiesUsed(template)) {
      grok.shell.info('Template is used in plates, cannot delete');
      return false;
    }
    await crud.deletePlateTemplate(template);
    return true;
  }
}
