import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import * as grok from "datagrok-api/grok";

export abstract class UaViewer {
  root: HTMLElement;
  name: string;
  setStyle: Function = null as any;
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

  protected constructor(name: string, setStyle?: Function | null, showName: boolean = true) {
    this.root = ui.div();
    this.name = name;

    if (setStyle)
      this.setStyle = setStyle;

    this.showName = showName;
  }

  reloadViewer() {
    this.root.innerHTML = '';
    let host = ui.block([]);
    if (this.setStyle)
      this.setStyle(host);
    if (this.showName)
      host.appendChild(ui.h1(this.name));
    let loader = ui.loader();
    host.appendChild(loader);

    this.setViewer(loader, host);

    this.root.append(host);
  }

  setViewer(loader: any, host: HTMLDivElement): void {
  }
}
