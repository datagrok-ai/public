import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {UaFilter} from "../filter2";
import {UaToolbox} from "../ua-toolbox";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {UaQueryViewer} from "../viewers/abstract/ua-query-viewer";

export abstract class UaView extends DG.ViewBase {
  uaToolbox: UaToolbox;
  viewers: UaQueryViewer[] = [];

  protected constructor(viewName: string, uaToolbox: UaToolbox) {
    super();
    this.name = viewName;
    this.uaToolbox = uaToolbox;
    this.toolbox = uaToolbox.rootAccordion.root;

    this.initViewers();
  }

  abstract initViewers() : any; // for sync and async code
}