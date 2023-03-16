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

  protected constructor(name: string, setStyle?: Function | null, showName: boolean = false) {
    this.root = ui.box();
    this.name = name;

    if (setStyle)
      this.setStyle = setStyle;

    this.showName = showName;
  }

  reloadViewer(header: boolean = true) {
    this.root.innerHTML = '';
    const host = ui.box();
    if (this.setStyle)
      this.setStyle(host);
    const loader = ui.loader();
    if (header) {
      const nameDiv = ui.divH([], {style: {alignItems: 'end'}});
      nameDiv.style.flexGrow = '0';
      if (this.showName)
        nameDiv.append(ui.h1(this.name, {style: {margin: '15px 0 3px 15px'}}));
      host.appendChild(nameDiv);
      this.setViewer(loader, host, nameDiv);
    } else this.setViewer(loader, host);
    host.appendChild(loader);
    this.root.append(host);
  }

  setViewer(loader: any, host: HTMLDivElement, nameDiv?: HTMLElement): void {}
}
