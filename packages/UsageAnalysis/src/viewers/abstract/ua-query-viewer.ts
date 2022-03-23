import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import * as grok from "datagrok-api/grok";
import {UaViewer} from "./ua-viewer";

export abstract class UaQueryViewer extends UaViewer{

  queryName: string;
  viewerFunction: Function;
  staticFilter: Object = {};
  filter: Object = {};

  protected constructor(name: string, queryName: string, viewerFunction: Function,
                        setStyle?: Function | null, staticFilter?: Object | null, filter?: UaFilter | null, showName: boolean = true) {
    super(name, setStyle, showName);

    this.queryName = queryName;
    this.viewerFunction = viewerFunction;

    if (staticFilter)
      this.staticFilter = staticFilter;
    if (filter)
      this.filter = filter;

    this.init();
  }

  setViewer(loader: any, host: HTMLDivElement) {
    let filter = {...this.filter, ...this.staticFilter}

    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      host.appendChild(this.viewerFunction(dataFrame));
      host.removeChild(loader);
    });
  }

  reload(filter: UaFilter) {
  };

  init() : void {
  }

}
