import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {toDart} from "datagrok-api/dg";
import * as rxjs from 'rxjs';
import * as GridUtils from '../utils/GridUtils';

const VIEWERS : Array<DG.Viewer> = [];

function indexofViewer(viewer: DG.Viewer) : number {
  for (let n = 0; n < VIEWERS.length; ++n) {
    if(toDart(VIEWERS[n]) === toDart(viewer))
      return n;
  }
  return -1;
}

let handlerTooltip : rxjs.Subscription | null = null

export class TooltipManager {

  private static init() : void {
    grok.events.onViewerClosed.subscribe((args: any) => {
      const viewer : DG.Viewer = (args as DG.ViewerArgs).viewer;
      const idx = indexofViewer(viewer);
      if(idx >= 0)
        VIEWERS.splice(idx, 1);

      if (VIEWERS.length === 0) {
        handlerTooltip?.unsubscribe();
        handlerTooltip = null;
      }
    });

    let x = 0;
    let y = 0;
    rxjs.fromEvent<MouseEvent>(document, 'mousemove').subscribe((e : MouseEvent) => {
      x = e.clientX;
      y = e.clientY;
    });

    handlerTooltip = grok.events.onTooltipShown.subscribe((args) => {
      //const context = args.args.context;
      const element : HTMLElement = args.args.element;
      if (grok.shell.tv === null || grok.shell.tv === undefined || grok.shell.tv.viewers === undefined ||
        element.children.length != 1)
        return;

      const observerResize = new ResizeObserver( (entries: ResizeObserverEntry[]) => {
        const rcTT = element.getBoundingClientRect();
        for (const entry of entries) {
          //console.log(`The resize was modified.` + rcTT + ' ' + entry.contentRect);

          const eUnderTT = document.elementFromPoint(rcTT.x, rcTT.y);
          const viewers = Array.from(grok.shell.tv.viewers);
          let viewer = null;
          let v = null;
          for (let n = 0; n < viewers.length; ++n) {
            v = viewers[n];
            let rr = eUnderTT;
            while (rr !== null) {
              if (v.root === rr) {
                viewer = v;
                break;
              }
              rr = rr.parentElement;
            }
            if (viewer !== null)
              break;
          }

          if (viewer !== null && viewer.type == DG.VIEWER.GRID) {
            const grid : DG.Grid = viewer as DG.Grid;
            const rcGrid = grid.root.getBoundingClientRect();
            const cellGrid = grid.hitTest(x - rcGrid.x, y - rcGrid.y);
            if (cellGrid === null || cellGrid === undefined)
              return;

            const renderer = GridUtils.getGridColumnRenderer(cellGrid.gridColumn);
            if(renderer === null)
              return;

            const tooltip = renderer.tooltip(cellGrid);
            if (tooltip === null)
              return;

            while (element.firstChild) {
              element.removeChild(element.firstChild);
            }

            element.appendChild(tooltip);
            args.preventDefault();
          }
        }
      });
      observerResize.observe(element);
    });
  }


  static isRegisted(viewer: DG.Grid) : boolean {
    const idx = indexofViewer(viewer);
    return idx >= 0;
  }

  static register(viewer: DG.Grid) : void {
    if (indexofViewer(viewer) >= 0)
      return;

    if (handlerTooltip === null) {
      TooltipManager.init()
    }

    VIEWERS.push((viewer));
  }
}
