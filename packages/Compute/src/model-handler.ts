import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from "wu";

export class ModelHandler extends DG.ObjectHandler {
  get type() {
    return 'Model';
  }

  async getById(id: string): Promise<DG.Func> {
    return await grok.dapi.functions.find(id);
  }

  async refresh(x: DG.Script): Promise<DG.Func> {
    return await this.getById(x.id);
  }

  getLanguageIcon(language: string) {
    if (language == 'grok')
      return ui.iconSvg('project');
    return ui.iconImage('script', `/images/entities/${language}.png`);
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof DG.Func && x.hasTag('model');
  }

  renderIcon(x: DG.Func, context: any = null): HTMLElement {
    if (x.options['icon'] != null && (x.options['icon'].startsWith('http://') || x.options['icon'].startsWith('https://'))) {
      return ui.iconImage('model-icon', x.options['icon']);
    }
    if (x instanceof DG.Script) {
      return this.getLanguageIcon(x.language);
    }
    let iconUrl = x.package.getIconUrl();
    if (x.options['icon'] != null) {
      let packagePathSegments = iconUrl.split('/');
      packagePathSegments.pop();
      packagePathSegments.push(x.options['icon']);
      iconUrl = packagePathSegments.join('/')
    }
    return ui.iconImage(x.package.name, iconUrl);
  }

  renderMarkup(x: DG.Func): HTMLElement {
    let markup = ui.span([this.renderIcon(x), ui.label(x.friendlyName)]);
    let c = grok.functions.getCurrentCall();
    markup.ondblclick = (e) => {
      ModelHandler.openModel(x, c); }
    return markup;
  }

  renderProperties(x: DG.Func) {
    const a = ui.accordion('ComputeModel');
    a.context = x;
    a.addTitle(ui.span([this.renderIcon(x), ui.label(x.friendlyName), ui.contextActions(x), ui.star(x.id)]));
    a.addPane('Details', () => {
      return this.renderDetails(x);
    }, true);
    a.end();
    return a.root;
  }

  renderDetails(x: DG.Func) {
    return ui.divV([ui.markdown(x.description), ui.span([ui.render(x.author), ' created ', ui.render(x.createdOn)])], {style: {lineHeight: '150%'}});
  }

  renderTooltip(x: DG.Func) {
    let h = this.renderMarkup(x);
    h.classList.add('d4-link-label');
    let card = ui.bind(x, ui.divV([
      ui.h2(h),
      ui.divV([ui.divText(x.description, 'ui-description')], {style: {lineHeight: '150%'}})
    ]));
    let c = grok.functions.getCurrentCall();
    card.ondblclick = (e) => { ModelHandler.openModel(x, c); }
    return card;
  }

  renderCard(x: DG.Func, context?: any): HTMLElement {
    let h = this.renderMarkup(x);
    h.classList.add('d4-link-label');
    let card = ui.bind(x, ui.h2(h));
    let c = grok.functions.getCurrentCall();
    card.ondblclick = (e) => { ModelHandler.openModel(x, c); }
    return card;
  }

  static openModel(x: DG.Func, parentCall?: DG.FuncCall) {
    let fc = x.prepare();
    fc.aux['showOnTaskBar'] = false;

    function findModelCatalogView() {
      return wu(grok.shell.views).find((v) => v.parentCall?.func.name == 'modelCatalog');
    }

    function startModel(modelsView: DG.View | undefined, parentCall: DG.FuncCall | undefined) {
      if (parentCall != null) {
        parentCall.aux['view'] = modelsView;
        fc.parentCall = parentCall;
      }
      fc.edit();
    }

    let modelsView = findModelCatalogView();
    if (modelsView == null) {
      grok.functions.call('Compute:ModelCatalog').then((_) => {
        modelsView = findModelCatalogView();
        startModel(modelsView, modelsView?.parentCall);
      });
    } else {
      startModel(modelsView, parentCall);
    }
  }

  init() {
  }
}
