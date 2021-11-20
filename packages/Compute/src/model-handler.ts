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

  getLanguageIcon(language: string) {
    if (language == 'grok')
      return ui.iconSvg('project');
    return ui.iconImage('script', `/images/entities/${language}.png`);
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof DG.Script && x.hasTag('model');
  }

  renderIcon(x: DG.Script, context: any = null): HTMLElement {
    return this.getLanguageIcon(x.language);
  }

  renderMarkup(x: DG.Script): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.name, {style: {marginLeft: '4px'}})], {style: {display: 'inline-flex'}});
  }

  renderProperties(x: DG.Script) {
    const a = ui.accordion();
    a.context = x;
    a.addTitle(ui.span([this.renderIcon(x), ui.label(x.friendlyName), ui.contextActions(x), ui.star(x.id)]));
    a.addPane('Details', () => {
      return this.renderDetails(x);
    });
    a.addCountPane('Usage', () => {
      return ui.span(['Usage statistics']);
    }, () => 99);
    a.end();
    return a.root;
  }

  renderDetails(x: DG.Script) {
    return ui.divV([ui.render(x.author), ui.render(x.createdOn)], {style: {lineHeight: '150%', marginTop: '16px'}});
  }

  renderTooltip(x: DG.Script) {
    return ui.divText(`${x.name}`);
  }

  renderCard(x: DG.Script, context?: any): HTMLElement {
    let card = ui.bind(x, ui.divV([
      ui.h2(this.renderMarkup(x)),
      ui.divText(x.description),
      this.renderDetails(x),
    ], 'd4-gallery-item'), {contextMenu: false});
    card.ondblclick = (e) => {
      this.openModel(x);
    }
    return card;
  }

  private openModel(x: DG.Script, parentView?: DG.ViewBase) {
    let pv = parentView ?? grok.shell.v;
    let view = DG.FunctionView.createFromFunc(x);
    if (pv.parentCall.func.name == 'modelCatalog' && pv instanceof DG.MultiView) {
      view.parentCall = pv.parentCall;
      view.parentView = pv;
      view.toolbox = ui.wait(async () => {
        let models = await grok.dapi.scripts
          .filter('#model')
          .list();
        return ui.divV(models.map((model) => ui.render(model, {onClick: (_) => this.openModel(model, pv)})), {style: {lineHeight: '150%'}});
      })
      //let mv: DG.MultiView = <DG.MultiView>grok.shell.v;
      //mv.addView(x.name, {factory: () => view, allowClose: true}, true);
    }
      grok.shell.addView(view);
  }

  init() {
    this.registerParamFunc('Open', async (t: DG.Script) => {

    });
  }
}
