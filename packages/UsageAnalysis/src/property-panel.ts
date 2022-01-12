import {UaFilterableViewer} from "./viewers/ua-filterable-viewer";
import {UaQueryViewer} from "./viewers/ua-query-viewer";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class PropertyPanel {
  private viewers: UaQueryViewer[];
  private acc: DG.Accordion;

  public constructor(viewers: UaQueryViewer[], accKey: string) {
    this.viewers = viewers;

    this.acc = DG.Accordion.create(accKey);

    for(let viewer of viewers) {
      this.acc.addPane(viewer.name, () => viewer.root);
    }
  }

  public getRoot() {
    return ui.block([this.acc.root]);
  }

}