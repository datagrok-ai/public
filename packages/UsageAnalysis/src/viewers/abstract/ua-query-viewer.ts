import * as ui from "datagrok-api/ui";
import {UaFilter} from "../../filter2";
import * as grok from "datagrok-api/grok";
import {UaViewer} from "./ua-viewer";
import * as DG from "datagrok-api/dg";

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

  setViewer(loader: any, host: HTMLDivElement, nameDiv: HTMLElement) {
    let filter = {...this.filter, ...this.staticFilter}

    grok.data.query('UsageAnalysis:' + this.queryName, filter).then((dataFrame) => {
      if (dataFrame.columns.byName('count') != null)
        dataFrame.columns.byName('count').tags['format'] = '#';

      let viewer: HTMLElement;
      if (dataFrame.rowCount > 0)
        viewer = this.viewerFunction(dataFrame);
      else
        viewer = ui.divText('No data', {style:{color:'var(--red-3)', paddingBottom:'25px'}})

      let grid = DG.Viewer.grid(dataFrame);
      grid.props.allowColSelection = false;
      let raw = false;

      let tableIcon = ui.iconFA('table', () => {
        if (!raw)
          $(viewer).replaceWith(grid.root);
        else
          $(grid.root).replaceWith(viewer);

        raw = !raw;
      });

      tableIcon.style.padding = '3px';
      tableIcon.style.margin = '0 3px';
      tableIcon.style.color = 'var(--grey-4)';

      tableIcon.addEventListener("mouseover",function(){tableIcon.style.color = 'var(--blue-1)'});
      tableIcon.addEventListener("mouseleave",function(){tableIcon.style.color = 'var(--grey-4)'});

      

      nameDiv.append(ui.tooltip.bind(tableIcon, 'Show grid'));

      host.appendChild(viewer);
      host.removeChild(loader);
    });
  }

  reload(filter: UaFilter) {
  };

  init() : void {
  }

}
