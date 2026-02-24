/* eslint-disable camelcase */
/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as rxjs from 'rxjs';
import * as PinnedUtils from '../pinned/PinnedUtils';
import {PinnedColumn} from '../pinned/PinnedColumn';
import {Grid, TableView} from 'datagrok-api/dg';

let m_dispatcher : MouseDispatcher | null = null;
function findPinnedColumn(eElem : Element) : PinnedColumn | null {
  if (!(eElem instanceof HTMLCanvasElement))
    return null;

  const itViews = grok.shell.views;
  const arViews = Array.from(itViews);
  let nPinnedColCount = -1;
  let view = null;
  let viewer = null;
  let itViewers = null;
  let arViewers = null;
  let colPinned = null;
  for (let nView=0; nView<arViews.length; ++nView) {
    view = arViews[nView];
    if (!(view instanceof TableView))
      continue;

    itViewers = view.viewers;
    arViewers = Array.from(itViewers);
    for (let nViewer=0; nViewer<arViewers.length; ++nViewer) {
      if (arViewers[nViewer].type != DG.VIEWER.GRID)
        continue;

      viewer = arViewers[nViewer];
      nPinnedColCount = PinnedUtils.getPinnedColumnCount(viewer as Grid);
      for (let nC=0; nC<nPinnedColCount; ++nC) {
        colPinned = PinnedUtils.getPinnedColumn(nC, viewer as Grid);
        if (colPinned.getRoot() === eElem)
          return colPinned;
      }
    }
  }

  return null;
}

function findElement(nXClient : number, nYClient : number) : Element | null {
  const eElem = document.elementFromPoint(nXClient, nYClient);

  return eElem;
}

export class MouseDispatcher {
  private m_eElementCurrent : Element | null;
  private m_colPinnedCurrent : PinnedColumn | null;
  private m_bDragging : boolean;

  private m_handlerMouseDown : rxjs.Subscription | null;
  private m_handlerMouseUp : rxjs.Subscription | null;
  private m_handlerContextMenu : rxjs.Subscription | null;
  private m_handlerMouseDblClk: rxjs.Subscription | null;
  private m_handlerMouseLeave : rxjs.Subscription | null;
  private m_handlerMouseMove : rxjs.Subscription | null;
  private m_handlerMouseWheel : rxjs.Subscription | null;
  private m_nXLastRightClick : number = -1;
  private m_nYLastRightClick : number = -1;

  constructor() {
    if (m_dispatcher != null)
      throw new Error('Mouse dispatcher can only have one instance.');

    m_dispatcher = this;

    this.m_eElementCurrent = null;
    this.m_colPinnedCurrent = null;

    this.m_bDragging = false;
    const dispatcherThis = this;

    this.m_handlerMouseDown = rxjs.fromEvent<MouseEvent>(document, 'mousedown').subscribe((e : MouseEvent) => {
      //we already know the element from the move
      dispatcherThis.m_bDragging = true;
      const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
      if (colPinnedCurr != null)
        colPinnedCurr.onMouseDown(e);
    });

    this.m_handlerMouseUp = rxjs.fromEvent<MouseEvent>(document, 'mouseup').subscribe((e : MouseEvent) => {
      dispatcherThis.m_bDragging = false;
      const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
      if (colPinnedCurr != null)
        colPinnedCurr.onMouseUp(e);
    });

    this.m_handlerContextMenu = rxjs.fromEvent<MouseEvent>(document, 'contextmenu').subscribe((e : MouseEvent) => {
      const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
      if (colPinnedCurr != null)
        colPinnedCurr.onContextMenu(e);
    });


    this.m_handlerMouseDblClk = rxjs.fromEvent<MouseEvent>(document, 'dblclick').subscribe((e : MouseEvent) => {
      //we already know the element from the move
      const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
      if (colPinnedCurr != null)
        colPinnedCurr.onMouseDblClick(e);
    });

    this.m_handlerMouseLeave = rxjs.fromEvent<MouseEvent>(document, 'mouseleave').subscribe((e : MouseEvent) => {
      const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
      if (colPinnedCurr != null)
        colPinnedCurr.onMouseLeave(e, false);

      dispatcherThis.setCurrentElement(null, null);
      //console.log('Mouse Left Document');
    });

    this.m_handlerMouseMove = rxjs.fromEvent<MouseEvent>(document, 'mousemove').subscribe((e : MouseEvent) => {
      if (dispatcherThis.m_bDragging) {
        const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
        if (colPinnedCurr != null)
          colPinnedCurr.onMouseDrag(e);

        return;
      }

      const eElem = findElement(e.clientX, e.clientY);
      //console.log(' Moving on Element ' + eElem);
      const eElemCurrent = dispatcherThis.getCurrentElement();
      if (eElem !== eElemCurrent) {
        const colPinnedCurr = dispatcherThis.getCurrentPinnedColumn();
        if (colPinnedCurr != null) {
          const rect = eElemCurrent == null ? null : eElemCurrent.getBoundingClientRect();
          const bOverlap = rect !== null && e.clientX >= rect.left && e.clientX <= rect.right && e.clientY >= rect.top && e.clientY <= rect.bottom;
          colPinnedCurr.onMouseLeave(e, bOverlap);
        }
        //console.log('Mouse Left Element ' + dispatcherThis.getCurrentElement());
        const colPinned = eElem === null ? null : findPinnedColumn(eElem);
        dispatcherThis.setCurrentElement(eElem, colPinned);
        if (colPinned !== null)
          colPinned.onMouseEnter(e);
        //console.log('Mouse on New Element ' + dispatcherThis.getCurrentElement() + ' ' + colPinned?.getGridColumn()?.name);
      }

      const colPinned = this.getCurrentPinnedColumn();
      if (colPinned !== null)
        colPinned.onMouseMove(e);
    });

    this.m_handlerMouseWheel = rxjs.fromEvent<WheelEvent>(document, 'wheel').subscribe((e : WheelEvent) => {
      const colPinned = this.getCurrentPinnedColumn();
      if (colPinned !== null) {
        colPinned.onMouseWheel(e);

        // event source is the document element
        //e.preventDefault();
        //e.stopPropagation();
        e.stopImmediatePropagation();
      }
    });
  }

  private getCurrentElement() : Element | null {
    return this.m_eElementCurrent;
  }

  private getCurrentPinnedColumn() : PinnedColumn | null {
    return this.m_colPinnedCurrent;
  }

  private setCurrentElement(eElem : Element | null, colPinned : PinnedColumn | null) {
    this.m_eElementCurrent = eElem;
    this.m_colPinnedCurrent = colPinned;
  }

  public dispose() : void {
    this.m_handlerMouseDown?.unsubscribe();
    this.m_handlerMouseDown = null;

    this.m_handlerMouseUp?.unsubscribe();
    this.m_handlerMouseUp = null;

    this.m_handlerContextMenu?.unsubscribe();
    this.m_handlerContextMenu = null;

    this.m_handlerMouseDblClk?.unsubscribe();
    this.m_handlerMouseDblClk = null;

    this.m_handlerMouseLeave?.unsubscribe();
    this.m_handlerMouseLeave = null;

    this.m_handlerMouseMove?.unsubscribe();
    this.m_handlerMouseMove = null;

    this.m_handlerMouseWheel?.unsubscribe();
    this.m_handlerMouseWheel = null;

    this.m_eElementCurrent = null;

    m_dispatcher = null;
  }

  static create() : MouseDispatcher {
    if (m_dispatcher != null)
      return m_dispatcher;

    return new MouseDispatcher();
  }
}
