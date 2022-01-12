import {UaFilterableViewer} from "./viewers/ua-filterable-viewer";
import {UaQueryViewer} from "./viewers/ua-query-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class PropertyPanel {
  private header: string;
  private viewers: UaQueryViewer[];
  private acc: DG.Accordion;

  public constructor(viewers: UaQueryViewer[], header: string, accordionKey: string) {
    this.viewers = viewers;
    this.header = header;

    this.acc = DG.Accordion.create(accordionKey);

    for(let viewer of viewers) {
      this.acc.addPane(viewer.name, () => viewer.root);
    }
  }

  public getRoot() {
    return  ui.block([
      ui.divV([
        ui.h1(this.header),
        this.acc.root
      ])
    ]);
  }

}