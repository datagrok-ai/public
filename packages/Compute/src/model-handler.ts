import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class ModelHandler extends DG.ObjectHandler {
  get type() {
    return 'Model';
  }

  async getById(id: string): Promise<DG.Script> {
    console.log('id:', id);
    return await grok.dapi.scripts.find(id);
  }

  async refresh(x: DG.Script): Promise<DG.Script> {
    let script = await this.getById(x.id);
    console.log('script:', script);
    return script;
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof DG.Script && x.hasTag('model');
  }

  renderIcon(x: DG.Script, context: any = null): HTMLElement {
    return ui.iconFA('function');
  }

  renderMarkup(x: DG.Script): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.name)]);
  }

  renderProperties(x: DG.Script) {
    const a = ui.accordion();
    a.addTitle(ui.span([this.renderIcon(x), ui.label(x.friendlyName), ui.contextActions(x), ui.star(x.id)]));
    a.addPane('Details', () => {
      return this.renderDetails(x);
    });
    a.addCountPane('Usage', () => {
      return ui.span(['Usage statistics']);
    }, () => 99);
    return a.root;
  }

  renderDetails(x: DG.Script) {
    return ui.tableFromMap({'Created': x.createdOn, 'By': x.author});
  }

  renderTooltip(x: DG.Script) {
    return ui.divText(`${x.name} is in the air!`);
  }

  renderCard(x: DG.Script, context?: any): HTMLElement {
    let card = ui.bind(x, ui.divV([
      ui.h2(x.friendlyName),
      this.renderDetails(x),
    ], 'd4-gallery-item'));
    card.ondblclick = (e) => {
      grok.shell.addView(DG.FunctionView.createFromFunc(x));
    }
    return card;
  }

  init() {
    this.registerParamFunc('Open', async (t: DG.Script) => {

    });
  }
}
