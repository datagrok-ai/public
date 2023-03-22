import * as ui from 'datagrok-api/ui';


export abstract class UaViewer {
  root: HTMLElement;
  name: string;
  setStyle: Function = null as any;

  static splineStyle: Object = {
    'aggrType': 'count',
    'innerChartMarginTop': 0,
    'innerChartMarginBottom': 0,
    'outerChartMarginTop': 5,
    'outerChartMarginBottom': 0,
    'yGlobalScale': false,
    'showTopPanel': false,
    'showMouseOverRowLine': false,
    'showXSelector': false,
    'showYSelectors': false,
    'showAggrSelectors': false,
    'showSplitSelector': false,
    // "showYAxis": false,
    'showMarkers': 'Never',
    'Title': 'Users',
  };

  static defaultBarchartOptions: Object = {
    valueAggrType: 'avg',
    style: 'dashboard',
  };

  static defaultChartOptions: Object = {
    style: 'dashboard',
  };

  protected constructor(name: string, setStyle?: Function | null) {
    this.root = ui.box();
    this.name = name;
    if (setStyle)
      this.setStyle = setStyle;
  }

  reloadViewer() {
    this.root.innerHTML = '';
    const host = ui.box();
    if (this.setStyle)
      this.setStyle(host);
    const loader = ui.loader();
    this.setViewer(loader, host);
    host.appendChild(loader);
    this.root.append(host);
  }

  setViewer(loader: any, host: HTMLDivElement): void {}
}
