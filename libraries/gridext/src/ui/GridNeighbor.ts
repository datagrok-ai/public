import * as DG from 'datagrok-api/dg';
import * as GridUtils from '../utils/GridUtils';
import * as ui from "datagrok-api/ui";
import * as grok from 'datagrok-api/grok';

export class GridNeighbor {

    private m_observerResizeGrid : ResizeObserver | null;
    private m_grid : DG.Grid | null;
    private m_root : HTMLElement | null;

    constructor(e : HTMLElement, grid : DG.Grid, nSize : number) {

        this.m_root  = e;
        this.m_grid = grid;
        const nW = nSize;
        const nHeight = grid.canvas.height;
        const tabIndex =  grid.canvas.getAttribute("tabIndex");
        if(tabIndex !== null)
        e.setAttribute("tabIndex", tabIndex);

        grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + "px";
        grid.overlay.style.left= (grid.overlay.offsetLeft + nW).toString() + "px";

        grid.canvas.style.width = (grid.canvas.offsetWidth - nW).toString() + "px";
        grid.overlay.style.width= (grid.overlay.offsetWidth - nW).toString() + "px";

        e.style.position = "absolute";
        e.style.left = 0 + "px";
        e.style.top = grid.canvas.offsetTop + "px";
        e.style.width = nW.toString() + "px";
        e.style.height = Math.round(nHeight/window.devicePixelRatio) + "px";

        if(grid.canvas.parentNode === null)
            throw new Error("Parent node for canvas cannot be null.");

        grid.canvas.parentNode.insertBefore(e, grid.canvas);

        const neighborThis = this;
        const eThis = e;
        this.m_observerResizeGrid = new ResizeObserver(function (entries : any) {
            const viewTable = grid.view;
            const bCurrent =  DG.toDart(grok.shell.v) === DG.toDart(viewTable);
            if(!bCurrent)
                return;

            if(grid.canvas.height !== eThis.offsetHeight) {
                eThis.style.top = grid.canvas.offsetTop + "px";
                eThis.style.width = nW + "px";
                eThis.style.height = Math.round(grid.canvas.height/window.devicePixelRatio) + "px";
            }
           neighborThis.onSizeChanged();
        });

        this.m_observerResizeGrid?.observe(grid.canvas);
    }

    get root() {return this.m_root};

    onSizeChanged() {}

    close() {

        if(this.m_grid === null || this.m_root === null)
            throw new Error('Grid or Root cannot be null.');

        this.m_observerResizeGrid?.disconnect();
        this.m_observerResizeGrid = null;

        this.m_grid.canvas.style.left = (this.m_grid.canvas.offsetLeft - this.m_root.offsetWidth).toString() + "px";
        this.m_grid.overlay.style.left= (this.m_grid.overlay.offsetLeft - this.m_root.offsetWidth).toString() + "px";
        this.m_grid.canvas.style.width = (this.m_grid.canvas.offsetWidth + this.m_root.offsetWidth).toString() + "px";
        this.m_grid.overlay.style.width= (this.m_grid.overlay.offsetWidth + this.m_root.offsetWidth).toString() + "px";

        if(this.m_root.parentNode !== null)
            this.m_root.parentNode.removeChild(this.m_root);

        this.m_root = null;
        this.m_grid = null;
    }
}