import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from "wu";
import {ModelCatalogView} from "./model-catalog-view";

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

  getVideoBackgroundUrl(x: DG.Func): string {
    let url = x.options['video-preview'];
    if (url == null)
      return '';
    if (url.startsWith('http://') || url.startsWith('https://'))
      return url;
    if (x instanceof DG.Script) {
      if (x.options['icon'] != null) {
        let iconPath = x.options['icon'].substring(0, x.options['icon'].lastIndexOf('/'));
        return `${iconPath}/${url}`;
      }
      return '';
    }
    let iconPath = x.package.getIconUrl();
    iconPath = iconPath.substring(0, iconPath.lastIndexOf('/'));
    return `${iconPath}/${url}`;
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
    let titleDiv = ui.div([
      ui.span([this.renderIcon(x), ui.label(x.friendlyName), ui.contextActions(x), ui.star(x.id)])]);
    a.addTitle(titleDiv);

    if (x.description != null)
      titleDiv.appendChild(ui.div([ui.markdown(x.description)], 'model-catalog-description'));
    titleDiv.appendChild(ui.tags(x));
    if ( x.options['dev.status'] != null)
      titleDiv.appendChild(ui.tableFromMap({Status: x.options['dev.status']}));

    if (x.options['video'] != null || x.options['education'] != null) {
      let educationLink = ui.link('Education materials', x.options['education']);
      ui.setDisplay(educationLink, x.options['education'] != null);
      let video = ui.div(ui.icons.play(() => {window.open(x.options['video'], '_blank')}), 'model-catalog-video');
      ui.setDisplay(video, x.options['video'] != null);
      video.style.backgroundImage = `url('${this.getVideoBackgroundUrl(x)}')`;

      a.addPane('Documentation', () => {
        return ui.divV([
          educationLink,
          video], {style: {lineHeight: '150%'}});
      }, true);
    }
    a.end();
    return a.root;
  }

  renderDetails(x: DG.Func) {
    let status = ui.markdown(`**Status:** ${x.options['dev.status']}`);
    ui.setDisplay(status, x.options['dev.status'] != null);
    return ui.divV([
      ui.markdown(`**Description: ** ${x.description}`),
      status,
      ui.markdown('**Tags:**'),
      ui.tags(x),
      ui.markdown('**Documentation:**'),
    ], {style: {lineHeight: '150%'}});
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

  static getModelCatalogView(): ModelCatalogView {
    let modelsView = this.findModelCatalogView();
    if (modelsView == null) {
      let mc: DG.Func = DG.Func.find({package:'Compute', name: 'ModelCatalog'})[0];
      let fc = mc.prepare();
      let view = new ModelCatalogView();
      view.name = 'Models';
      view.parentCall = fc;
      grok.shell.addView(view);
    }
    modelsView = this.findModelCatalogView();
    return modelsView as ModelCatalogView;
  }

  static findModelCatalogView(): ModelCatalogView | undefined {
    return wu(grok.shell.views).find((v) => v.parentCall?.func.name == 'modelCatalog') as ModelCatalogView;
  }

  static openModel(x: DG.Func, parentCall?: DG.FuncCall) {
    let fc = x.prepare();
    fc.edit();
  }

  static bindModel(fc: DG.FuncCall) {
    let modelsView = ModelHandler.getModelCatalogView();
    fc.aux['showOnTaskBar'] = false;
    if (modelsView != null) {
      let parentCall = modelsView.parentCall;
      if (parentCall != null) {
        parentCall.aux['view'] = modelsView;
        fc.parentCall = parentCall;
      }
    }
  }

  init() {
    setTimeout(() => {
      grok.functions.onBeforeRunAction.subscribe((fc) => {
        if (fc.func.hasTag('model') || fc.func.hasTag('model-editor')) {
          ModelHandler.bindModel(fc);
        }
      });
    });
  }
}
