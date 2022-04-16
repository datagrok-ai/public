import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelCatalogView} from "./model-catalog-view";
import {TYPE} from "datagrok-api/dg";
import {FunctionView} from '@datagrok-libraries/utils/src/function-view';

export class ModelHandler extends DG.ObjectHandler {
  get type() {
    return 'Model';
  }

  async getById(id: string): Promise<DG.Func> {
    console.log('id:', id);
    return await grok.dapi.functions.find(id);
  }

  async refresh(x: DG.Script): Promise<DG.Func> {
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
    return x instanceof DG.Func && x.hasTag('model');
  }

  renderIcon(x: DG.Func, context: any = null): HTMLElement {
    return x instanceof DG.Script ? this.getLanguageIcon(x.language) : ui.iconFA('lightning');
  }

  renderMarkup(x: DG.Func): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.friendlyName, {style: {marginLeft: '4px'}})], {style: {display: 'inline-flex'}});
  }

  renderProperties(x: DG.Func) {
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

  renderDetails(x: DG.Func) {
    return ui.divV([ui.render(x.author), ui.render(x.createdOn)], {style: {lineHeight: '150%'}});
  }

  renderTooltip(x: DG.Func) {
    return ui.divText(`${x.friendlyName}`);
  }

  renderCard(x: DG.Func, context?: any): HTMLElement {
    let card = ui.bind(x, ui.divV([
      ui.h2(this.renderMarkup(x)),
      ui.divText(x.description, 'ui-description'),
      this.renderDetails(x),
    ], 'd4-gallery-item'), {contextMenu: false});
    let c = grok.functions.getCurrentCall();
    card.ondblclick = (e) => {
      ModelHandler.openModel(x, c);
    }
    return card;
  }

  static openModel(x: DG.Func, parentCall?: DG.FuncCall) {
    let fc = x.prepare();
    if (parentCall != null)
      fc.parentCall = parentCall;
    fc.edit();
  }

  init() {
  }
}
