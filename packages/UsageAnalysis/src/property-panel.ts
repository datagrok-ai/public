import {UaFilterableQueryViewer} from "./viewers/ua-filterable-query-viewer";
import {UaQueryViewer} from "./viewers/abstract/ua-query-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {UaViewer} from "./viewers/abstract/ua-viewer";

export class PropertyPanel {
  private acc: DG.Accordion;
  private root: HTMLElement;

  public constructor(entity: DG.Entity | null, viewers: UaViewer[], header: string, accordionKey: string) {
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