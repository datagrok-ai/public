import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class CDMViewBase extends DG.ViewBase {

  constructor(name) {
    super({});
    this.name = name;
  }

  loaded = false;
  filterChanged = false;

  load(): void {
    this.createView();
    this.loaded = true;
  }

  createView() { }

  updateGlobalFilter() {}

  async propertyPanel() {
    return await this.createAccWithTitle(this.name).root;
  }

  createAccWithTitle(panelName: string, title = null) {
    let acc = ui.accordion(`${panelName} panel`);
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(`${title ?? panelName}`)]));
    return acc;
  }

}