import * as ui from "datagrok-api/ui";
import {UaFilter} from "../filter2";
import * as grok from "datagrok-api/grok";

export abstract class UaQueryViewer {
  root: HTMLElement;
  name: string;
  queryName: string;
  viewerFunction: Function;
  setStyle: Function = null as any;
  staticFilter: Object = {};
  filter: Object = {};
  showName: boolean;

  static splineStyle: Object = {
    "aggrType": "count",
    "innerChartMarginTop": 0,
    "innerChartMarginBottom": 0,
    "outerChartMarginTop": 5,
    "outerChartMarginBottom": 0,
    "yGlobalScale": false,
    "showTopPanel": false,
    "showMouseOverRowLine": false,
    "showXSelector": false,
    "showYSelectors": false,
    "showAggrSelectors": false,
    "showSplitSelector": false,
    "showYAxis": false,
    "showMarkers": "Never",
    "Title":"Users"
  };

  static defaultBarchartOptions: Object = {
    valueAggrType: 'avg',
    style: 'dashboard'
  };

  static defaultChartOptions: Object = {
    style: 'dashboard'
  };

  protected constructor(name: string, queryName: string, viewerFunction: Function,
                        setStyle?: Function, staticFilter?: Object, filter?: UaFilter, showName: boolean = true) {
    this.root = ui.div();
    this.name = name;
    this.queryName = queryName;
    this.viewerFunction = viewerFunction;

    if (setStyle)
      this.setStyle = setStyle;
    if (staticFilter)
      this.staticFilter = staticFilter;
    if (filter)
      this.filter = filter;

    this.showName = showName;

    this.init();
  }

  // addCardUsingDataframe(cardName: string, dataFrame: DG.DataFrame, viewer:any, supportUsers = true) {
  //     let host = ui.block([],'d4-item-card card');
  //     host.appendChild(ui.h1(cardName));
  //
  //     // if (cardName === 'Errors')
  //     //     grok.data.detectSemanticTypes(dataFrame);
  //     host.appendChild(viewer(dataFrame));
  //     return host;
  // }

  reloadViewer() {
    this.root.innerHTML = '';
    let host = ui.block([]);
    if (this.setStyle)
      this.setStyle(host);
    if (this.showName)
      host.appendChild(ui.h1(this.name));
    let loader = ui.loader();
    host.appendChild(loader);

    let filter = {...this.filter, ...this.staticFilter}

    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      // if (cardName === 'Errors')
      //     grok.data.detectSemanticTypes(dataFrame);
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';
      host.appendChild(this.viewerFunction(dataFrame));
      host.removeChild(loader);
    });
    this.root.append(host);
  }

  init() : void {
  }

  reload(filter: UaFilter) {
  };
}
