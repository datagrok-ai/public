import * as DG from "datagrok-api/dg";
import {ErrorUtils} from "../../utils/ErrorUtils";
import {DefaultColumnCtx} from "./ColumnCtx";
import {IncompatibleFilterError} from "./IncompatibleFilterError";
import {DGUtils} from "../../utils/DGUtils";

export class FiltersUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}


FiltersUtils.FILTER_CSS_CLASS_NAME = "d4-filter";
FiltersUtils.FILTER_LABEL_CSS_CLASS_NAME = "d4-filter-column-name";
FiltersUtils.FILTER_PANEL_CSS_CLASS_NAME = "d4-root d4-filter-group d4-viewer d4-filters ui-box";
FiltersUtils.FILTER_CSS_CLASS_NAMEE =  "single-col-query-device";


FiltersUtils.addFilter = function(filter, viewerFilters) {
    ErrorUtils.verifyClass(filter, DG.Filter);

    let eFilterGroup = viewerFilters.root.querySelector('.d4-filter-group-wrapper');
    if(eFilterGroup === null)
     eFilterGroup = viewerFilters.root;

    const root = filter.root;
    const lstClassesCss = root.classList;

    lstClassesCss.remove("ui-div");
    lstClassesCss.add("d4-flex-col");
    lstClassesCss.add("d4-filter");

    const lstChildren = eFilterGroup.children;
    if(lstChildren.length === 0)
     eFilterGroup.appendChild(filter.root);
    else
     eFilterGroup.insertBefore(filter.root, lstChildren[0]);

    setTimeout(() => {
        filter.root.style = "";
    }, 500);
}


FiltersUtils.removeFilter = function(filter, viewerFilters) {
    ErrorUtils.verifyClass(filter, DG.Filter);

    let eFilterGroup = viewerFilters.root.querySelector('.d4-filter-group-wrapper');
    if(eFilterGroup === null)
        eFilterGroup = viewerFilters.root;

    const lstChildren = eFilterGroup.children;

    let eChild = null;
    for (let i = 0; i < lstChildren.length; i++) {
        eChild = lstChildren[i];
        if(eChild === filter.root)
        {
            eFilterGroup.removeChild(eChild);
            return true;
        }
    }
    return false;
}




FiltersUtils.getFiltersViewer = function()
{
    const viewCurrent = grok.shell.v;
    const viewTable = viewCurrent;//grok.shell.getTableView(dframe.name);

    if(viewCurrent.dart === viewTable.dart)
    {
        let ccc = 0;
    }

    const ar = DGUtils.findViewersbyType(viewTable, DG.VIEWER.FILTERS);
    if(ar.length > 0)
     return ar[0];
    /*
    const lstViewers = viewTable.viewers;
    for(var viewer of lstViewers)
    {
        if(viewer.type === DG.VIEWER.FILTERS)
            return viewer;
    }     */

    return null;
}


FiltersUtils.hasCompatibleFilter = function(col)
{
    return col.dart.m_bHasFilter;
}


FiltersUtils.setHasCompatibleFilter = function(col, bHasFilter)
{
    col.dart.m_bHasFilter = bHasFilter;
}

FiltersUtils.createFilter = function(col, arFilterClasses)
{
    ErrorUtils.verifyClass(col, DG.Column);

    const ctxCol = new DefaultColumnCtx(col);
    let filter = null;

    const nDeviceCount = arFilterClasses.length;
    for (let nDevice = 0; nDevice < nDeviceCount; ++nDevice) {
        try {filter = new arFilterClasses[nDevice](ctxCol);}
        catch(e)
        {
          if(e instanceof IncompatibleFilterError)
            continue;
          else throw new Error(e);
        }
        break;
    }

    return filter;
}





FiltersUtils.getFiltersPanel = function(dframe)
{
    if(dframe === undefined)
        throw new Error("Data frame cannot be undefined");

    let viewerFilter = FiltersUtils.getFiltersViewer();

    if(viewerFilter === null)
        viewerFilter = FiltersUtils.openFilterViewer(dframe);

    if(viewerFilter === null)
        return null;

   return viewerFilter.dart.m_panelFilters;
}


 /*
FiltersUtils.isFilterPanelOpened = function()
{
    const c = document.getElementsByClassName(FiltersUtils.FILTER_PANEL_CSS_CLASS_NAME);
    return  c.length > 0;
}  */

FiltersUtils.m_fnFilterOpener = null;
FiltersUtils.m_fnFilterCloser = null;

FiltersUtils.getFilterOpener = function() {
    return FiltersUtils.m_fnFilterOpener;
}

FiltersUtils.setFilterOpener = function(fnOpener) {
    FiltersUtils.m_fnFilterOpener = fnOpener;
}

FiltersUtils.getFilterCloser = function() {
    return FiltersUtils.m_fnFilterCloser;
}

FiltersUtils.setFilterCloser = function(fnCloser) {
    FiltersUtils.m_fnFilterCloser = fnCloser;
}


FiltersUtils.openFilterViewer = function(dframe)
{
    //if(dframe === undefined)
      //  throw new Error("Data frame cannot be undefined");
    ErrorUtils.verifyClass(dframe, DG.DataFrame);

    let viewerFilters = null;

    if(FiltersUtils.m_fnFilterOpener !== null)
        viewerFilters = FiltersUtils.m_fnFilterOpener();
    else {
        const viewSS = grok.shell.v;//my changes getTableView(dframe.name);
        viewerFilters = DG.Viewer.filters(dframe);
        viewSS.addViewer(viewerFilters);
    }
    /*
    try {grok.shell.dockManager.dock(viewerFilters.root, 'left', null, 'Query Devices');}
    catch (e)
    {
      return null;
    } */


    return viewerFilters;
}




FiltersUtils.openFilter = function(column) {

    const dframe = column.dataFrame;
    const panelFilters = FiltersUtils.getFiltersPanel(dframe);
    if (panelFilters === null)
        return null;

        let filter = panelFilters.getFilter(column.name);
        if (filter === null)
            filter = panelFilters.addFilter(column);

        return filter;

    //return null;
}



FiltersUtils.getFilter = function(renderer, colGrid)
{
    const column = renderer.adjustColumn(colGrid.column);
    const dframe = column.dataFrame;
    const panelFilters = FiltersUtils.getFiltersPanel(dframe);
    let filter = null;
    if(panelFilters !== null) {
        const colGrid = cell.gridColumn;
        filter = panelFilters.getFilter(column.name);
        if (filter === null)
            filter = panelFilters.addFilter(column);
    }
    return filter;
}



FiltersUtils.floatFilterToColumn = function(renderer, cell)
{
    const column = renderer.adjustColumn(cell);
    const dframe = column.dataFrame;
    const panelFilters = FiltersUtils.getFiltersPanel(dframe);
    if(panelFilters !== null)
    {
        const colGrid = cell.gridColumn;
        let filter = panelFilters.getFilter(column.name);
        if (filter === null)
            filter = panelFilters.addFilter(column);

        if(filter === null)
        {
            alert("There is no compatible filter for column '" + cell.tableColumn.name + "'");
            return;
        }

        const rc = cell.bounds;
        const nH = cell.grid.getOptions().look.colHeaderHeight;
        rc.height = nH;

        const rcGrid = cell.grid.canvas.getBoundingClientRect();
        filter.float( rcGrid.left + rc.x, rcGrid.top + rc.height);

        const nColGrid = colGrid.idx;
    }
}



   /*
FiltersUtils.getFilterPanelRoot = function(strColName)
{
    let c = document.getElementsByClassName(FiltersUtils.FILTER_PANEL_CSS_CLASS_NAME);
    if(c.length > 0)
    {
        const eDivPanel = c.item(0);
        return eDivPanel;
    }

    c = document.getElementsByClassName(FiltersUtils.FILTER_CSS_CLASS_NAME);
    if(c.length > 0)
    {
        const eDivPanel = c.item(0);
        return eDivPanel;
    }


    return null;
}    */


FiltersUtils.getFilterRootbyColName = function(strColName)
{
    const c = document.getElementsByClassName(FiltersUtils.FILTER_LABEL_CSS_CLASS_NAME);

    let e = null;
    for(var n=0; n<c.length; ++n)
    {
        e = c.item(n);
        if(e.innerHTML === strColName)
        {
            const eP = e.parentElement.parentElement;
            return eP;
        }
        // e.style.display = "none";
    }

    return null;
}


FiltersUtils.showFilters = function(bShow) {
    setTimeout(function() {
        const c = document.getElementsByClassName(FiltersUtils.FILTER_CSS_CLASS_NAME);
        let e = null;
        for (var n = 0; n < c.length; ++n) {
            e = c.item(n);
            e.style.display = bShow ? "block" : "none";
        }
    },500);
}
