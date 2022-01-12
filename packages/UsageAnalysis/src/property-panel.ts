import {UaFilterableViewer} from "./viewers/ua-filterable-viewer";
import {UaQueryViewer} from "./viewers/ua-query-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class PropertyPanel {
  private acc: DG.Accordion;
  private root: HTMLElement;

  public constructor(entity: DG.Entity | null, viewers: UaQueryViewer[], header: string, accordionKey: string) {
    this.acc = DG.Accordion.create(accordionKey);

    if (entity != null)
      this.acc.addPane('Details', () => ui.render(entity));

    for(let viewer of viewers) {
      this.acc.addPane(viewer.name, () => viewer.root);
    }

    this.root = ui.block([
      ui.divV([
        ui.h1(header),
        this.acc.root
      ])
    ]);
  }

  public getRoot(): HTMLElement {
    return this.root;
  }

}