// Defines the way Datagrok handles entities of the specified type
import * as DG from "datagrok-api/dg";
import {Template} from "./template";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

export class TemplateHandler extends DG.ObjectHandler {
  get type() {
    return 'Template'
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof Template;
  }

  /*
    getCanvasRenderer(x) { return new FruitCanvasRenderer(); }
    getGridCellRenderer(x) { return new FruitGridCellRenderer(); }*/

  renderIcon(x: Template, context: any = null): HTMLElement {
    return ui.iconFA('apple-alt');
  }

  renderMarkup(x: Template): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.name)]);
  }

  renderProperties(x: Template) {
    return ui.divText(`Properties for ${x.name}`);
  }

  renderTooltip(x: Template) {
    return ui.divText(`${x.name} is in the air!`);
  }

  renderCard(x: Template, context?: any): HTMLElement {
    return ui.bind(x, ui.divV([
      this.renderMarkup(x),
      ui.divText(`Context: ${context}`)
    ], 'd4-gallery-item'));
  }

  init() {
    this.registerParamFunc('Open', async (t: Template) => {
      t.context = DG.Context.create();
      t.context.setVariable('template', t);
      let query = await grok.functions.eval(t.query);
      let call = query.prepare(t.queryParams);
      call.context = t.context;

      console.log ('handler', call.context.getVariable('template'));
      await call.call(true, undefined, {processed: false});
      t.data = call.getOutputParamValue();
     // console.log(t.data);
    //  grok.shell.addTableView(t.data!);
    });
  }
}