import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Observable, Subject} from 'rxjs';

import {ILogger} from '../utils/logger';

export class GridNeighbor {
  private m_observerResizeGrid: ResizeObserver | null;
  private m_grid: DG.Grid | null;
  private m_root: HTMLElement | null;

  constructor(e: HTMLElement, grid: DG.Grid, nSize: number,
    private readonly logger?: ILogger, resizeable: boolean = true
  ) {
    const logPrefix = `${this.toLog()}.ctor()`;
    this.logger?.debug(`${logPrefix}, start`);

    this.m_root = e;
    this.m_grid = grid;
    const nW = nSize;
    const nHeight = grid.canvas.height;
    const tabIndex = grid.canvas.getAttribute('tabIndex');
    if (tabIndex !== null)
      e.setAttribute('tabIndex', tabIndex);

    grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + 'px';
    grid.overlay.style.left = (grid.overlay.offsetLeft + nW).toString() + 'px';

    grid.canvas.style.width = (grid.canvas.offsetWidth - nW).toString() + 'px';
    grid.overlay.style.width = (grid.overlay.offsetWidth - nW).toString() + 'px';

    e.style.position = 'absolute';
    e.style.left = 0 + 'px';
    e.style.top = grid.canvas.offsetTop + 'px';
    e.style.width = nW.toString() + 'px';
    e.style.height = Math.round(nHeight / window.devicePixelRatio) + 'px';

    if (grid.canvas.parentNode === null)
      throw new Error('Parent node for canvas cannot be null.');

    grid.canvas.parentNode.insertBefore(e, grid.canvas);

    const neighborThis = this;
    const eThis = e;

    function adjustHorzScroll() {
      const horzScroll = grid.root.querySelector('.d4-grid-horz-scroll') as HTMLElement;
      if (horzScroll) {
        const neighborWidth = eThis.offsetWidth ?? 0;
        horzScroll.style.left = `${neighborWidth}px`;
      }
    }
    this.m_observerResizeGrid = new ResizeObserver(function(entries: any) {
      requestAnimationFrame(() => {
        adjustHorzScroll();
        const viewTable = grid.view;
        const bCurrent = DG.toDart(grok.shell.v) === DG.toDart(viewTable);
        if (!bCurrent)
          return;

        if (grid.canvas.height !== Math.floor(eThis.offsetHeight * window.devicePixelRatio) ) {
          eThis.style.top = grid.canvas.offsetTop + 'px';
          // eThis.style.width = nW + 'px';
          eThis.style.height = Math.round(grid.canvas.height / window.devicePixelRatio) + 'px';
        }
        neighborThis.onSizeChanged();
      });
    });
    adjustHorzScroll();
    this.m_observerResizeGrid?.observe(grid.canvas);
    if (resizeable)
      attachResizer(e, grid, this);
  }

  private static objCounter: number = -1;
  private readonly objId: number = ++GridNeighbor.objCounter;

  private toLog(): string { return `GridNeighbor<${this.objId}>`; }

  get root() { return this.m_root; };

  private m_onClosed: Subject<void> = new Subject<void>();

  get onClosed(): Observable<void> { return this.m_onClosed; }

  onSizeChanged() {}

  close() {
    const logPrefix = `${this.toLog()}.close()`;
    this.logger?.debug(`${logPrefix}, start`);

    if (this.m_grid === null || this.m_root === null)
      throw new Error('Grid or Root cannot be null.');

    this.m_observerResizeGrid?.disconnect();
    this.m_observerResizeGrid = null;

    this.m_grid.canvas.style.left = (this.m_grid.canvas.offsetLeft - this.m_root.offsetWidth).toString() + 'px';
    this.m_grid.overlay.style.left = (this.m_grid.overlay.offsetLeft - this.m_root.offsetWidth).toString() + 'px';
    this.m_grid.canvas.style.width = (this.m_grid.canvas.offsetWidth + this.m_root.offsetWidth).toString() + 'px';
    this.m_grid.overlay.style.width = (this.m_grid.overlay.offsetWidth + this.m_root.offsetWidth).toString() + 'px';

    if (this.m_root.parentNode !== null)
      this.m_root.parentNode.removeChild(this.m_root);

    this.m_onClosed.next();

    // reset horz scroll
    const horzScroll = this.m_grid?.root.querySelector('.d4-grid-horz-scroll') as HTMLElement;
    if (horzScroll)
      horzScroll.style.removeProperty('left');

    this.m_root = null;
    this.m_grid = null;

    this.logger?.debug(`${logPrefix}, end`);
  }
}


function attachResizer(e: HTMLElement, grid: DG.Grid, gn: GridNeighbor) {
  let dragging = false;
  let dragStartX = 0;

  e.addEventListener('mousemove', (ev) => {
    // e.style.cursor = 'col-resize';
    const rect = e.getBoundingClientRect();
    const rectLeft = rect.left;
    const mousePos = ev.clientX - rectLeft;

    // console.log(`mousePos: ${mousePos}, e.offsetWidth: ${e.offsetWidth}`);
    if (mousePos > e.offsetWidth - 5)
      e.style.cursor = 'col-resize';
    else
      e.style.removeProperty('cursor');
  });
  e.addEventListener('mouseleave', (ev) => {
    if (!dragging)
      e.style.removeProperty('cursor');
  });

  e.addEventListener('mousedown', (ev) => {
    const rect = e.getBoundingClientRect();
    const rectLeft = rect.left;
    const mousePos = ev.clientX - rectLeft;
    if (mousePos > e.offsetWidth - 5) {
      dragging = true;
      dragStartX = ev.clientX;
      e.style.cursor = 'col-resize';
      ev.stopImmediatePropagation();
      ev.preventDefault();
    }
  });

  function onMouseUp(_ev: MouseEvent) {
    dragging = false;
  }

  function onMouseMove(ev: MouseEvent) {
    if (!dragging)
      return;
    const dx = ev.clientX - dragStartX;
    dragStartX = ev.clientX;
    const newWidth = e.offsetWidth + dx;
    if (newWidth < 100 || newWidth > 800)
      return;
    e.style.width = newWidth.toString() + 'px';
    grid.canvas.style.left = (grid.canvas.offsetLeft + dx).toString() + 'px';
    grid.overlay.style.left = (grid.overlay.offsetLeft + dx).toString() + 'px';
    grid.canvas.style.width = (grid.canvas.offsetWidth - dx).toString() + 'px';
    grid.overlay.style.width = (grid.overlay.offsetWidth - dx).toString() + 'px';
    gn.onSizeChanged();
    grid.invalidate();
  }

  window.addEventListener('mouseup', onMouseUp);
  window.addEventListener('mousemove', onMouseMove);
  const onCloseSub = gn.onClosed.subscribe(() => {
    window.removeEventListener('mouseup', onMouseUp);
    window.removeEventListener('mousemove', onMouseMove);
    onCloseSub.unsubscribe();
  });
}
