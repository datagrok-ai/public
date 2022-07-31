// Defines the way Datagrok handles entities of the specified type
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {HitTriageTemplate} from "./hit-triage-app";

export class TemplateHandler extends DG.ObjectHandler {
  get type() {
    return 'HitTriageTemplate'
  }

  // Checks whether this is the handler for [x]
  isApplicable(x: any) {
    return x instanceof HitTriageTemplate;
  }

  renderIcon(x: HitTriageTemplate, context: any = null): HTMLElement {
    return ui.iconFA('apple-alt');
  }

  renderMarkup(x: HitTriageTemplate): HTMLElement {
    return ui.span([this.renderIcon(x), ui.label(x.project)]);
  }

  renderProperties(x: HitTriageTemplate) {
    return ui.divText(`Properties for ${x.project}`);
  }

  renderTooltip(x: HitTriageTemplate) {
    return ui.divText(`${x.project} is in the air!`);
  }

  renderCard(x: HitTriageTemplate, context?: any): HTMLElement {
    return ui.bind(x, ui.divV([
      this.renderMarkup(x),
      ui.divText(`Context: ${context}`)
    ], 'd4-gallery-item'));
  }

  init() {
    this.registerParamFunc('Open', async (t: HitTriageTemplate) => {

    });
  }
}