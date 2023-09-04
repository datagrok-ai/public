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

  renderIcon(x: DG.Func, context: any = null): HTMLElement {
    if (x.options['icon'] != null && (x.options['icon'].startsWith('http://') || x.options['icon'].startsWith('https://'))) {
      return ui.iconImage('model-icon', x.options['icon']);
    }
    if (x instanceof DG.Script) {
      return this.getLanguageIcon(x.language);
    }
    x = DG.Func.find({package: x.package.name, name: x.name})[0];
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
    markup.ondblclick = (e) => { ModelHandler.openModel(x); }
    markup.onclick = (e) => { ModelHandler.openHelp(x); }
    return markup;
  }

  async renderPreview(x: DG.Func) {
    const editorName = x.options.editor ?? 'Compute:RichFunctionViewEditor';
    const editor = await grok.functions.find(editorName);
    if (editor !== null && editor instanceof DG.Func) {
      const viewCall = editor.prepare({'call': x.prepare()});
      await viewCall.call(false, undefined, {processed: true});
      const view = viewCall.getOutputParamValue();
      if (view instanceof DG.View)
        return view;
    }
    //@ts-ignore
    return super.renderPreview(x);
  }

  renderProperties(func: DG.Func) {
    const a = ui.accordion('ComputeModel');
    a.context = func;
    let titleDiv = ui.div([
      ui.span([this.renderIcon(func), ui.label(func.friendlyName), ui.contextActions(func), ui.star(func.id)])]);
    a.addTitle(titleDiv);

    if (func.description != null)
      titleDiv.appendChild(ui.div([ui.markdown(func.description)], 'model-catalog-description'));
    titleDiv.appendChild(ui.tags(func));
    if ( func.options['dev.status'] != null)
      titleDiv.appendChild(ui.tableFromMap({Status: func.options['dev.status']}));

    a.end();
    return a.root;
  }

  renderTooltip(x: DG.Func) {
    let h = this.renderMarkup(x);
    h.classList.add('d4-link-label');
    let card = ui.bind(x, ui.divV([
      ui.h2(h),
      ui.divV([ui.divText(x.description, 'ui-description')], {style: {lineHeight: '150%'}})
    ]));
    card.ondblclick = (e) => { ModelHandler.openModel(x); }
    return card;
  }

  renderCard(x: DG.Func, context?: any): HTMLElement {
    let h = this.renderMarkup(x);
    h.classList.add('d4-link-label');
    let card = ui.bind(x, ui.h2(h));
    card.ondblclick = (e) => { ModelHandler.openModel(x); }
    card.onclick = (e) => { ModelHandler.openHelp(x); }
    return card;
  }

  static async openHelp(func: DG.Func) {
    if (func.options['readme'] != null) {
      const path = `System:AppData/${func.package.name}/${func.options['readme']}`;
      const readmeText = await grok.dapi.files.readAsText(path);

      grok.shell.windows.help.showHelp(ui.markdown(readmeText));
      grok.shell.windows.showHelp = true;
    }
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

  static openModel(x: DG.Func) {
    let fc = x.prepare();
    fc.edit();
  }

  static openModelFromFunccall(funcCall: DG.FuncCall) {
    funcCall.edit();
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
