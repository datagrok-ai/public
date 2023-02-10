import * as ui from 'datagrok-api/ui';

export abstract class UaViewer {
  root: HTMLElement;
  name: string;
  setStyle: Function = null as any;
  showName: boolean;

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

  protected constructor(name: string, setStyle?: Function | null, showName: boolean = true) {
    this.root = ui.div();
    this.name = name;

    if (setStyle)
      this.setStyle = setStyle;

    this.showName = showName;
  }

  reloadViewer() {
    this.root.innerHTML = '';
    const host = ui.block([]);
    if (this.setStyle)
      this.setStyle(host);

    const nameDiv = ui.divH([], {style: {alignItems: 'center'}});
    if (this.showName)
      nameDiv.append(ui.h1(this.name, {style: {margin: '6px 0'}}));
    host.appendChild(nameDiv);

    const loader = ui.loader();
    host.appendChild(loader);

    this.setViewer(loader, host, nameDiv);

    this.root.append(host);
  }

  setViewer(loader: any, host: HTMLDivElement, nameDiv: HTMLElement): void {
  }
}
