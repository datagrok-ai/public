import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

let handlerViewerAddedDock = null;
let handlerViewLayoutApplyingDock = null;
let handlerViewLayoutAppliedDock = null;
let layoutManager = null;

/**
 * The TableViewLayoutManager class defines a top level interface for layout manager
 * for DG viewers that places them in a TableView in accordance with a specific
 * dockng pattern.
 */
export class TableViewLayoutManager {
    /**
     * Docks a specified viewer in the current TableView.
     * This method is called by the framework. Users shouldn't call it directly.
     * To be implemented in sub-classes.
     * @param viewer the viewer to be docked.
     */
    layout(viewer) {}
    //to be implemented in sub-classes
}

/**
 * Uninstalls the currently installed layout manager in the DG platform.
 * This method should be called before applying any saved TableView layouts.
 * @returns {boolean} true if the manager was uninstalled successfully, otherwise false.
 */
TableViewLayoutManager.uninstall = function () {
    return TableViewLayoutManager.install(null);
};

function _uninstall() {
    if (handlerViewerAddedDock !== null) {
        handlerViewerAddedDock.unsubscribe();
        handlerViewerAddedDock = null;

        //layoutManager = null;
        return true;
    }

    return false;
}


function _install() {
    if (layoutManager === null) {
        return false;
    }
    handlerViewerAddedDock = grok.events.onViewerAdded.subscribe((args) => {
        const { viewer } = args.args;

        layoutManager.layout(viewer);
    });

    return true;
}

/**
 * Installs a specified layout manager into the DG platform.
 * This method should be called after any saved TableView layouts have been applied.
 * @param layout the specified layout manager to be installed, or null to uninstalled the
 * currently installed layout.
 * @returns {boolean} true if the manager was installed successfully, otherwise false.
 */
TableViewLayoutManager.install = function (layout) {
    if (layout === null) {
        if (handlerViewerAddedDock !== null) {
            handlerViewerAddedDock.unsubscribe();
            handlerViewerAddedDock = null;
        }

        if (handlerViewLayoutApplyingDock !== null) {
            handlerViewLayoutApplyingDock.unsubscribe();
            handlerViewLayoutApplyingDock = null;
        }

        if (handlerViewLayoutAppliedDock !== null) {
            handlerViewLayoutAppliedDock.unsubscribe();
            handlerViewLayoutAppliedDock = null;
        }

        layoutManager = null;
        return true;
    }

    if (!(layout instanceof TableViewLayoutManager)) { return false; }

    layoutManager = layout;
    handlerViewerAddedDock = grok.events.onViewerAdded.subscribe((args) => {
        const { viewer } = args.args;
        layoutManager.layout(viewer);
    });


    if (handlerViewLayoutApplyingDock === null) {
        handlerViewLayoutApplyingDock = grok.events.onViewLayoutApplying.subscribe((layout) => {
            _uninstall();
        });
    }

    if (handlerViewLayoutAppliedDock === null) {
        handlerViewLayoutAppliedDock = grok.events.onViewLayoutApplied.subscribe((layout) => {
            _install();
        });
    }

    return true;
};

/**
 * Finds and returns the rightmost viewers in the TableView.
 * @param ar a returned array that after the cal will contain the rightmost viewers.
 * @param itViewers an iterator that contains viewers from the TableView.
 * @param viewerExclude the specified viewer to be excluded from the search. This is the viewer to be docked.
 * @param grid the grid viewer.
 */
function fillRightmostViewers(ar, itViewers, viewerExclude, grid) {
    while (ar.length > 0) {
        ar.pop();
    }
    /*
  const arViewers = Array.from(itViewers);

  const v1 = grok.shell.v;
  const v2 = grok.shell.v;

  const b = v1 === v2;
  const bb = DG.toDart(v1) === DG.toDart(v2);
      */
    const rcGrid = grid.root.getBoundingClientRect();
    const nXGrid = rcGrid.x + rcGrid.width - 1;

    let nX = -1;
    let nXRight = -1;
    let rc = null;

    const it = itViewers[Symbol.iterator]();
    let ob = it.next();
    let v = null;
    while (!ob.done) {
        v = ob.value;
        if (v.type !== DG.VIEWER.FILTERS && v.type !== DG.VIEWER.GRID && DG.toDart(viewerExclude) !== DG.toDart(v)
            && v.root.offsetParent !== null) {
            rc = v.root.getBoundingClientRect();
            nX = rc.x + rc.width - 1;
            if (nX > nXGrid) {
                if (nXRight < 0 || nX === nXRight) {
                    ar.push(v);
                    nXRight = nX;
                } else if (nX > nXRight) {
                    while (ar.length > 0) {
                        ar.pop();
                    }
                    ar.push(v);
                    nXRight = nX;
                }
            }
        }
        ob = it.next();
    }
}

/**
 * The FlowTileLayoutManager provides implementation for layout manager based on
 * a docking algorithm that places viewers on right side of a TableView using
 * vertical tiling pattern of a specified vertical size.
 */
export class FlowTileLayoutManager extends TableViewLayoutManager {
    /**
     * Constructs a new instance of the FlowTileLayoutManager class for a specified vertical size.
     * @param nRowCount the maximum number of tiles in a vertical column.
     */
    constructor(nRowCount) {
        super();

        this.m_nRowCount = nRowCount === undefined ? 2 : nRowCount;
        this.m_fRatioHorz = 0.11;
        this.m_fRatioVert = 0.5;
    }

    /**
     * Returns the maximum number of tiles in a vertical column.
     * @returns {number} an integer representing the maximum number of tiles.
     */
    getRowCount() { return this.m_nRowCount; }

    /**
     * Returns the horizontal ratio that determines the width of the docked viewers
     * relative to the horizontal size of the TableView.
     * Please see the details in the GROK API https://datagrok.ai/js-api/classes/dg.dockmanager#dock
     * @returns {number} a floating point number in the range [0,1] representing the ratio.
     */
    getHorzRatio() { return this.m_fRatioHorz; }

    /**
     * Sets the horizontal ratio that determines the width of the docked viewers
     * relative to the horizontal size of the TableView.
     * Please see the details in the GROK API https://datagrok.ai/js-api/classes/dg.dockmanager#dock
     * @returns {number} a floating point number in the range [0,1] representing the ratio.
     */
    setHorzRatio(fRatio) {
        this.m_fRatioHorz = fRatio;
    }

    /**
     * Returns the vertical ratio that determines the height of the docked viewers
     * relative to the vertical size of the TableView.
     * Please see the details in the GROK API https://datagrok.ai/js-api/classes/dg.dockmanager#dock
     * @returns {number} a floating point number in the range [0,1] representing the ratio.
     */
    getVertRatio() { return this.m_fRatioVert; }

    /**
     * Sets the vertical ratio that determines the height of the docked viewers
     * relative to the vertical size of the TableView.
     * Please see the details in the GROK API https://datagrok.ai/js-api/classes/dg.dockmanager#dock
     * @returns {number} a floating point number in the range [0,1] representing the ratio.
     */
    setVertRatio(fRatio) {
        this.m_fRatioVert = fRatio;
    }


    /**
     * Overrides the implementation of the parent class.
     * @param viewer the viewer to be docked.
     */
    layout(viewer) {
        const { type } = viewer;
        if (type === DG.VIEWER.GRID) {
            return;
        }

        const { grid } = viewer.view;

        const manager = viewer.view.dockManager;

        const nodeViewer = manager.findNode(viewer.root);
        manager.close(viewer.root);
        nodeViewer.container.destroy();

        if (type === DG.VIEWER.FILTERS) {
            manager.dock(viewer.root, DG.DOCK_TYPE.LEFT, manager.rootNode, type, this.m_fRatioHorz);
            return;
        }

        const viewTable = viewer.view;
        const itViewers = viewTable.viewers;

        const ar = [];
        fillRightmostViewers(ar, itViewers, viewer, grid);

        if (ar.length === 0 || ar.length >= this.m_nRowCount) {
            manager.dock(viewer.root, DG.DOCK_TYPE.RIGHT, manager.rootNode, type, this.m_fRatioHorz);
        } else {
            const node = manager.findNode(ar[ar.length - 1].root);
            manager.dock(viewer.root, DG.DOCK_TYPE.DOWN, node, type, this.m_fRatioVert);
        }
    }
}
