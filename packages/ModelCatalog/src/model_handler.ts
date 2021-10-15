import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";

export class ModelHandler extends DG.ObjectHandler {
  get type() {
    return 'Model'
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof DG.Script && x.hasTag("modelhub");
  }

  renderIcon(x: DG.Script, context: any = null): HTMLElement {
    return ui.iconFA('function');
  }

  renderMarkup(x: DG.Script): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.name)]);
  }

  renderProperties(x: DG.Script) {
    let a = ui.accordion();
    a.addTitle(ui.span([this.renderIcon(x), ui.label(x.friendlyName), ui.contextActions(x), ui.star(x.id)]));
    a.addPane('Details', () => {return this.renderDetails(x)});
    a.addCountPane('Usage', () => {return ui.span(['Usage statistics'])}, () => 99);
    return a.root;
  }

  renderDetails(x: DG.Script) {
    return ui.tableFromMap({'Created': x.createdOn, 'By': x.author});
  }

  renderTooltip(x: DG.Script) {
    return ui.divText(`${x.name} is in the air!`);
  }

  renderCard(x: DG.Script, context?: any): HTMLElement {
    return ui.bind(x, ui.divV([
      ui.h2(x.friendlyName),
      this.renderDetails(x)
    ], 'd4-gallery-item'));
  }

  init() {
    this.registerParamFunc('Open', async (t: DG.Script) => {

    });
  }
}