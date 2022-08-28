import * as DG from "datagrok-api/dg";
import {SemType} from "../entity/SemType";
import {PrimitiveSemType} from "../entity/primitive/PrimitiveSemType";
import {PrimitiveQNumSemType} from "../entity/primitive/PrimitiveQNumSemType";
import {PrimitiveDatetSemType} from "../entity/primitive/PrimitiveDatetSemType";
import {PrimitiveStringSemType} from "../entity/primitive/PrimitiveStringSemType";
import {FiltersUtils} from "../ui/filters/FiltersUtils";
import {render} from "datagrok-api/ui";
import {ErrorUtils} from "./ErrorUtils";
import {CpdSemType} from "../entity/gnf/cpd/CpdSemType";
/**
 * The GridUtils class contains advanced functionality for DataGrok Grid Viewer.
 * This class is intended for use by the application framework and contains
 * experimental functionality that relies on the DataGrok new feature development and
 * might not exist in future versions. Please use this class on your own risk.
 */
export class GridUtils {
    constructor() {
        throw new Error("Never create instances of this class");
    }
}

GridUtils.isTooltipEnabled = function (grid)
{
    return grid.dart.m_bTooltipEnabled === undefined ? false : grid.dart.m_bTooltipEnabled;
}

GridUtils.enableTooltip = function (grid, bEnable)
{
 grid.dart.m_bTooltipEnabled = bEnable;
}

GridUtils.getRowsHeader = function(grid, strColName)
{
    if(grid.dart.m_mapRowHeaders === undefined)
        return null;

    const header = grid.dart.m_mapRowHeaders.get(strColName);
    if(header === null || header === undefined)
        return null;

    return header;
}

GridUtils.getRowsHeaderCount = function(grid)
{
    if(grid.dart.m_mapRowHeaders === undefined)
        return 0;

    return grid.dart.m_mapRowHeaders.size;
}

GridUtils.setRowsHeader = function(grid, header, bClose)
{
    if(grid.dart.m_mapRowHeaders === undefined)
       grid.dart.m_mapRowHeaders = new Map();

    if(bClose) {
        grid.dart.m_mapRowHeaders.delete(header.getColumnName());
        return;
    }

    grid.dart.m_mapRowHeaders.set(header.getColumnName(), header);
}


GridUtils.repaintRowsHeaders = function(grid)
{
    if(grid.dart.m_mapRowHeaders === undefined ||  grid.dart.m_mapRowHeaders === null)
        return;

    const it = grid.dart.m_mapRowHeaders.values();
    for(let header of it)
    {
        if(header === null || header === undefined)
            continue;

        header.invalidate();
    }
}


GridUtils.setRenderersConfig = function(grid, config)
{
    ErrorUtils.verifyClass(grid, DG.Grid);
    grid.dart.m_configRenderers = config;
}

GridUtils.getRenderersConfig = function(grid)
{
    ErrorUtils.verifyClass(grid, DG.Grid);
    if(grid.dart.m_configRenderers === undefined)
        return null;

    return grid.dart.m_configRenderers;
}



GridUtils.getSortOptionContext = function(grid)
{
    return grid.dart.m_ctxSort;
}

GridUtils.setSortOptionContext = function(grid, ctx)
{
    grid.dart.m_ctxSort = ctx;
}

GridUtils.getCellBackgroundColor = function(cellGrid)
{
    ErrorUtils.verifyClass(cellGrid, DG.GridCell);
    const cell = cellGrid.cell;
    ErrorUtils.verifyClass(cell, DG.Cell);

    let cr = null;
    try {
        cr = DG.Color.getCellColorHtml(cell);
    } catch (e) {
        console.error("ERROR: Couldn't get cell color html");
        cr = null;
    }

    return cr;
}

GridUtils.repaint = function(grid)
{
    const scrollH = grid.horzScroll;
    const nMin = scrollH.min;
    const nMax = scrollH.max;

    scrollH.setValues(scrollH.minRange, scrollH.maxRange, nMin+1, nMax);
    scrollH.setValues(scrollH.minRange, scrollH.maxRange, nMin, nMax);

    //setInterval(() => { grid.invalidate(); }, 1000);
}


GridUtils.sort = async function (cell, nRecord, nOption)
{
    const colGrid = cell.gridColumn;
    const grid = cell.grid;
    const col = colGrid.column;
    const typeSem = GridUtils.getSemType(col);
    if(nOption >= typeSem.getSortOptionCount())
        throw new Error("Sorting option (" + nOption + ") is out of bounds [0, " + (typeSem.getSortOptionCount() -1) + "]");

    const nRecordCount = col.length;
    const arData = new Array(nRecordCount);
    const ctxSort = GridUtils.getSortOptionContext(grid);

    let bAscend = true;
    if(typeSem.supportsSortDirectionSwitch(nOption) && ctxSort !== null && ctxSort.m_nIdxColGrid === colGrid.idx &&  ctxSort.m_nSortOption === nOption)
        bAscend = !ctxSort.m_bSortAscend;

    await SemType.sort(arData, col, nRecord, nOption, bAscend);
    GridUtils.setRowOrder(grid, arData);
    GridUtils.setSortOptionContext(grid, new SortOptionCtx(colGrid.idx, nOption, bAscend));
}



GridUtils.setRowOrder = function(grid, arIndices) {

    const scrollH = grid.horzScroll;
    const nMin = scrollH.min;
    const nMax = scrollH.max;

    grid.setRowOrder(arIndices);
    scrollH.setValues(scrollH.minRange, scrollH.maxRange, nMin, nMax);
    GridUtils.repaintRowsHeaders(grid);
}


GridUtils.assignSemType = function(col)
{
    let typeSem = GridUtils.getSemType(col);
    if(typeSem !== null && typeSem !== undefined)
        return typeSem;

    const type = col.type;
    if(type === DG.COLUMN_TYPE.QNUM)
        typeSem = PrimitiveQNumSemType.Instance;
    else if(type === DG.COLUMN_TYPE.DATE_TIME)
        typeSem = PrimitiveDatetSemType.Instance;
    else if(type === DG.COLUMN_TYPE.STRING)
        typeSem = PrimitiveStringSemType.Instance;
    else
        typeSem = PrimitiveSemType.Instance;

    GridUtils.setSemType(col, typeSem);

    return typeSem;
}


GridUtils.assignSemTypes = function(dframe)
{
    let col = null;
    let type = null;
    let typeSem = null;
    const nColCount = dframe.columns.length;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        col = dframe.columns.byIndex(nCol);
        typeSem = GridUtils.assignSemType(col);//GridUtils.getSemType(col);
            /*
        if(typeSem !== null && typeSem !== undefined)
            continue;

        type = col.type;
        if(type === DG.TYPE.QNUM)
          typeSem = PrimitiveQNumSemType.Instance;
        else if(type === DG.TYPE.DATE)
            typeSem = PrimitiveDatetSemType.Instance;
        else if(type === DG.TYPE.STRING)
            typeSem = PrimitiveStringSemType.Instance;
        else
            typeSem = PrimitiveSemType.Instance;

        GridUtils.setSemType(col, typeSem);*/
    }
}



GridUtils.index4ColName = function(arCols, strName) {

    let col = null;
    for (var nC = 0; nC < arCols.length; ++nC) {
        if(arCols[nC].name === strName)
         return nC;
    }

    return -1;
}

GridUtils.names2Cols = function(dframe, arColNames) {

    ErrorUtils.verifyNotUndefinedNotNull(dframe);
    ErrorUtils.verifyNotUndefinedNotNull(arColNames);

    let col = null;
    const arCols = new Array(arColNames.length);
    for (var nC = 0; nC < arColNames.length; ++nC) {
        col = dframe.columns.byName(arColNames[nC]);
        arCols[nC] = col;
    }

    return arCols;
}

GridUtils.col2Names = function(arCols) {

    ErrorUtils.verifyNotUndefinedNotNull(arCols);

    let col = null;
    const arColNames = new Array(arCols.length);
    for (var nC = 0; nC < arCols.length; ++nC) {
        col = arCols[nC];
        arColNames[nC] = col.name;
    }

    return arColNames;
}


GridUtils.getSemType = function(col) {
    const typeSem = col.dart.m_typeSem;
    return typeSem === undefined ? null : typeSem;
}


GridUtils.setSemType = function(col, typeSem) {
    col.dart.m_typeSem = typeSem;
}



GridUtils.findColumnBySemType = function(dframe, typeSem) {
    let column = null;
    let nColCount = dframe.columns.length;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        column = dframe.columns.byIndex(nCol);
        if(GridUtils.getSemType(column) instanceof typeSem)
            return nCol;
    }

    return -1;
}

GridUtils.findGridColumnBySemType = function(grid, typeSem) {
    let cc = grid.constructor.name;

    let arCols = grid.columns.dart.list;
    if(arCols === undefined)
        arCols = grid.columns.dart.b;

    let column = null;
    let nColCount = arCols.length;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        column = grid.columns.byIndex(nCol);
        if(column.column !== null && GridUtils.getSemType(column.column) instanceof typeSem)
         return nCol;
    }

    return -1;
}

GridUtils.findPrevGridColumn = function(nColGrid, grid, bVisible)
{

    const lstCols = grid.columns;
    let col = null;
    for(let nCol=nColGrid -1; nCol>=0; --nCol)
    {
        col = lstCols.byIndex(nCol);
        if((col.visible && bVisible) || (!col.visible && !bVisible))
            return col;
    }

    return null;
}



GridUtils.calcTextrHeight = function(g, strText, font, nWidth, nInsets) {

    g.font = font;
    const tm = g.measureText("W");
    const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nInsets;

    const arLines = [];
    GridUtils.calcWordsCellLayout(g, arLines, strText, nWidth);
    const nHeight = arLines.length*nHFont;

    return nHeight + nInsets;
}




//Evaluates the preferred height for a Grid's table header.
//ctx the canvas' graphics context. Should be initialized with font and alignment flags before calling this function.
//df the DataGrok data frame
//nWidth the specifide column width
//nInsets the specified gap between the text lines as well as between the cell's border and text. Usually 1-2 pixels.
//return the height value as integer value in pixels
GridUtils.calcPrefColHeaderHeight = function(ctx, df, nWidth, nInsets) {

    const tm = ctx.measureText("W");
    //const nYInset = 2;
    const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nInsets;

    let col = null;
    const lstCols = df.columns.toList();
    const nColCount = lstCols.length;

    const arLines = [];   //reused instance
    let nHeight = 0;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        col = lstCols[nCol];
        GridUtils.calcWordsCellLayout(ctx, arLines, col.name, nWidth);
        nHeight = Math.max(nHeight, arLines.length*nHFont);
    }

    return nHeight + nInsets;
}


GridUtils.calcPrefColHeaderHeightNew = function(g, dframe, configGridRenderers, nWidthDefault, nInsets) {

    const nYInsets = 2; //distance between lines
    g.font = "13px Arial";
    const tm = g.measureText("W");
    const nYInset = 2;
    const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nInsets;

    let col = null;
    const lstCols = dframe.columns.toList();
    const nColCount = lstCols.length;

    let i = null;
    let nWidth = -1;
    let typeSem = null;
    let renderer = null;
    const arLines = [];   //reused instance
    let nHeight = 0;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        col = lstCols[nCol];
        typeSem = GridUtils.getSemType(col);
        if(typeSem === null || typeSem === undefined)
         renderer = configGridRenderers.getDefaultColHeaderRenderer();
        else
         renderer = configGridRenderers.getColHeaderRenderer(typeSem.constructor);

        if(renderer !== null)
        {
            i = renderer.getInsets();
            nWidth = renderer.getPreferredCellWidth() - (i.getL() + i.getR());
        }
        else nWidth = nWidthDefault;

        GridUtils.calcWordsCellLayout(g, arLines, col.name, nWidth);
        nHeight = Math.max(nHeight, arLines.length*(nHFont+nYInsets)+ (renderer !== null ? i.getT() + i.getB() : 0));
    }

    return nHeight + nInsets;
}



//Calculates individual words layout from an arbitrary string within a specified bounded range.
//The text is split only by whitespace as other characters are considered to a part of individual words.
//The priority is given to keep the words as a whole. Truncation occurs only used when an individual word cannot fit withing the range's bounds
//ctx the canvas' graphics context. Should be initialized with font and alignment flags before calling this function.
//arLines an array that after the call will contain the words layout. Each element of the array is another array representing the layout of an individual line.
//strText a string containing text to be layouted.
//nWidth the width of bounded range.
GridUtils.calcWordsCellLayout = function(ctx, arLines, strText, nWidth) {

    while(arLines.length > 0) {
        arLines.pop();
    }

    if(strText === "")
        return;

    arLines.push([]);

    const arWords = [];
    GridUtils.splitWithWhitespaces(strText, arWords);

    let tm = null;
    let strWord = "";
    let nLine = 0;
    let nWLine = 0;
    let nWWord = -1;
    let nCharCount = -1;
    for(let nWord=0; nWord<arWords.length; ++nWord) {
        strWord = arWords[nWord];
        tm = ctx.measureText(strWord);
        nWWord = tm.width;
        if(nWLine + nWWord > nWidth)  //the word cannot fit within the cell's bounds
        {
            if(arLines[nLine].length === 0)  // the word gets truncated if the line starts with it. Otherwise it goes to the next line (eventually it will be truncated)
            {
                nCharCount = GridUtils.fitWordPart(ctx, strWord, nWidth);
                if(nCharCount === 0)
                    return;

                strWord = strWord.substring(0, nCharCount);
                arLines[nLine].push(strWord);   //replace the word with the remaining part to repeat on the next iteration

                arWords[nWord] = arWords[nWord].substring(nCharCount);
            }

            --nWord;
            nWLine = 0;
            arLines.push([]);
            ++nLine;
            continue;
        }

        arLines[nLine].push(strWord);
        nWLine += nWWord;
    }
}

//Split an arbitrary text into words by whitespace character. Whitespaces as considered as individual words, and thus are included into the output.
//strText a string containing some text to be split.
//ar the specified array that after the call will contain individual words including whitespaces as legitimate words.
GridUtils.splitWithWhitespaces = function(strText, ar) {

    while(ar.length > 0) {
        ar.pop();
    }

    var nIdxStart=0;
    var nIdxSpace=-1;

    while(nIdxStart < strText.length)
    {
        nIdxSpace = strText.indexOf(" ", nIdxStart);
        if(nIdxSpace >= 0)
        {
            if(nIdxSpace> nIdxStart)
                ar.push(strText.substring(nIdxStart, nIdxSpace))

            ar.push(" ");
            nIdxStart = nIdxSpace + 1;
        }
        else
        {
            ar.push(strText.substring(nIdxStart, strText.length));
            break;
        }
    }
}


//Calculates the number of characters from a word that fits within a specified bounded range.
//strWord the specified word to be split
//nWidth the width of the specified bounded range.
GridUtils.fitWordPart = function(g, strWord, nWidth) {
    let str = "";
    let nWWWord = -1;
    let nLength = strWord.length;
    for (var n = nLength; n >= 0; --n) {
        str = strWord.substring(0, n);
        nWWWord = g.measureText(str).width;
        if(nWWWord <= nWidth)
            return n;
    }

    return 0;
}


GridUtils.adjustRowSizeToColumnWidth = function(grid) {

    const nColCpd = GridUtils.findGridColumnBySemType(grid, CpdSemType);
    const colGrid = grid.columns.byIndex(nColCpd); //row header
    const nH = Math.floor(colGrid.width/grid.dart.m_fZoomValue);
    GridUtils.setPreferredRowHeight(grid, nH);
    GridUtils.setRowHeight(grid, colGrid.width);
    /*
    grid.setOptions({
        rowHeight: colGrid.width
    });*/

    //colGrid = grid.columns.byIndex(0);
    //colGrid.width = 150;
}

GridUtils.initZoom = function(grid)
{
    grid.dart.MIN_ZOOM_VALUE = 0.35;
    grid.dart.MAX_ZOOM_VALUE = 4.0;
    grid.dart.getMinimumZoomValue = function() {return grid.dart.MIN_ZOOM_VALUE;}
    grid.dart.getMaximumZoomValue = function() {return grid.dart.MAX_ZOOM_VALUE;}

    grid.dart.DEFAULT_PREFERRED_ROW_HEIGHT = 25;
    grid.dart.DEFAULT_PREFERRED_COLUMN_WIDTH = 75;

    grid.dart.MIN_CELL_SIZE = 15;
    grid.dart.MIN_COLUMN_WIDTH = 0;

    grid.dart.m_fZoomValue = NaN;
    grid.dart.getZoomValue = function()
    {
        let sss = 0;
        return this.m_fZoomValue;
    }

    grid.dart.setZoomValue = function(fValue)
    {
        if(fValue < this.MIN_ZOOM_VALUE || fValue > this.MAX_ZOOM_VALUE)
            throw new Eoor("Zoom value is out of bounds [" + this.MIN_ZOOM_VALUE + "," +
                this.MAX_ZOOM_VALUE +
                "] " + fValue);


        let nValueOld = this.m_fZoomValue;
        if(nValueOld === fValue)
            return;

        this.m_fZoomValue = fValue;
        let fK = fValue;

        let nHeight = grid.dart.m_nHMinimum + Math.floor((grid.dart.m_nHPreferred - grid.dart.m_nHMinimum)*fK);
        if(nHeight < this.MIN_CELL_SIZE)
            nHeight = this.MIN_CELL_SIZE;

        if(nHeight !== GridUtils.getRowHeight(grid))
        {
            GridUtils.setRowHeight(grid, nHeight);
         /*
            grid.setOptions({
                rowHeight: nHeight
            });*/

        }

        let nColCount = GridUtils.getColumnCount(grid);
        let colGrid = null;
        let colGridPrev = null;
        let nWidth = -1;
        let bChangeCols = false;
        let nX = -1;

        for(var nCol=1; nCol<nColCount; ++nCol)
        {
            colGrid = grid.columns.byIndex(nCol);
            colGridPrev = GridUtils.findPrevGridColumn(nCol, grid, true);//grid.columns.byIndex(nCol-1);
            //if(colGridPrev == null) {
                //colGridPrev = grid.columns.byIndex(nCol-2);
            //}

            nWidth = colGrid.dart.m_nWMinimum + Math.floor((colGrid.dart.m_nWPreferred - colGrid.dart.m_nWMinimum)*fK);

            if(nWidth < this.MIN_CELL_SIZE)
                nWidth = this.MIN_CELL_SIZE;

            if(nWidth !== colGrid.dart.m_nWPreferred)
            {
                GridUtils.setColumnWidth(colGrid, nWidth);
                nX = colGridPrev === null ? 0 : GridUtils.getColumnX(colGridPrev) + GridUtils.getColumnWidth(colGridPrev);
                GridUtils.setColumnX(colGrid, nX);
                bChangeCols = true;
            }
        }

        if(bChangeCols) {
            //Repaint has been done alredy when setting column width
            //colGrid = grid.columns.byIndex(0); //header
            //let nWTmp =  GridUtils.getColumnWidth(colGrid);//colGrid.d.width;
            //GridUtils.setColumnWidth(colGrid, nWTmp -1);//colGrid.d.width = nWTmp-1;
            //colGrid.width = nWTmp;
        }
    }
}





GridUtils.getZoomValue = function(grid)
{
    return grid.dart.getZoomValue();
}

GridUtils.setZoomValue = function(grid, fZoomValue)
{
    grid.dart.setZoomValue(fZoomValue);
}

GridUtils.getRowHeight = function(grid) {
    const options = grid.getOptions();
    const nH =  options.look.rowHeight;
    if(nH === undefined)//DG bug
        return 0;

    return nH;
}


GridUtils.setRowHeight = function(grid, nH) {

    grid.setOptions({
        rowHeight: nH
    });

    GridUtils.repaintRowsHeaders(grid);
}


GridUtils.getPreferredRowHeight = function(grid)
{
    return grid.dart.m_nHPreferred;
}

GridUtils.setPreferredRowHeight = function(grid, nH)
{
    grid.dart.m_nHPreferred = nH;
}

GridUtils.getMinimumRowHeight = function(grid)
{
    return grid.dart.m_nHMinimum;
}

GridUtils.setMinimumRowHeight = function(grid, nH)
{
    grid.dart.m_nHMinimum = nH;
}

GridUtils.getColumnCount = function(grid) {
    let arCols = grid.columns.dart.list;
    if (arCols === undefined)
        arCols = grid.columns.dart.b;

    const nColCount = arCols.length;
    return nColCount;
}


GridUtils.getPreferredColumnWidth = function(colGrid)
{
    return colGrid.dart.m_nWPreferred;
}

GridUtils.setPreferredColumnWidth = function(colGrid, nW)
{
    colGrid.dart.m_nWPreferred = nW;
}

GridUtils.getMinimumColumnWidth = function(colGrid)
{
    return colGrid.dart.m_nWMinimum;
}

GridUtils.setMinimumColumnWidth = function(colGrid, nW)
{
    colGrid.dart.m_nWMinimum = nW;
}


GridUtils.getColumnHeaderY = function(grid)
{
    let nY = 0;//my changes grid.d._colLabelsBox.top;
    return nY;
}

GridUtils.getColumnHeaderHeight = function(grid)
{
    const options = grid.getOptions(true);
    const nHColHeader = options.look.colHeaderHeight;
    return nHColHeader;
}

 /*
GridUtils.getRowHeight = function(grid)
{
    const options = grid.getOptions();
    const nHRow = options.look.rowHeight;
    return nHRow;
}  */

GridUtils.getColumnX = function(colGrid)
{
    if(colGrid.dart.left !== undefined)
        return colGrid.dart.left;

    const nX =  colGrid.left;             //null bug in DG, temporary
    const nW = colGrid.width;
    //if(nX === undefined)
    //  nX = colGrid.d.fx;
    return nX;
}

GridUtils.setColumnX = function(colGrid, nX)
{
    colGrid.dart.left = nX;
    //colGrid.left = nX;
}


GridUtils.getColumnWidth = function(colGrid)
{
    //if(colGrid.d.width === undefined)
    //  return colGrid.d.fr;

    return colGrid.width;
}

GridUtils.setColumnWidth = function(colGrid, nW)
{
    /*
    if(colGrid.d.width === undefined)
        colGrid.d.fr = nW;
    else
        colGrid.d.width = nW;
      */
    colGrid.width = nW;
}


GridUtils.isFirstViewportColumn = function(col, grid)
{
 const scroll = grid.horzScroll;
 const nMin = scroll.min;
 const b = col.left <= nMin && nMin <= col.left + col.width;
 return b;
}


GridUtils.fillVisibleGridCells = function(arColRowIdxs, grid)
{
    if(arColRowIdxs.length !== 4)
        throw new Error("Array to cobtain bound row column indices must have the length 4.");

    const arRowIdxs = [];
    const arColIdxs = [];
    const lstVisibleCells = grid.getVisibleCells();
    for(var cellGTmp of lstVisibleCells)
    {
        if(!arRowIdxs.includes(cellGTmp.gridRow))
         arRowIdxs.push(cellGTmp.gridRow);

        if(!arColIdxs.includes(cellGTmp.gridColumn.idx))
         arColIdxs.push(cellGTmp.gridColumn.idx);
    }
       /*
    const cellHeaderTopLeft = grid.hitTest(3, 3);
    const colTopLeft = cellHeaderTopLeft.gridColumn;
    const strNameTL = colTopLeft.name;

    const scrollH = grid.horzScroll;
    const minR = scrollH.minRange;
    const maxR = scrollH.maxRange;
    const min = scrollH.min;
    const max = scrollH.max;

    let col = null;
    let nWCols = colTopLeft.width - (scrollH.min - colTopLeft.left);
    const nColCount = grid.columns.length;
    let nColTopRight = colTopLeft.idx;
    for(nColTopRight=colTopLeft.idx+1; nColTopRight<nColCount; ++nColTopRight)
    {
        if(nWCols >= grid.canvas.width)
        {
            --nColTopRight;
            //console.log("top right " + col.name);
            break;
        }

        col = grid.columns.byIndex(nColTopRight);
        if(!col.visible)
            continue;

        nWCols += col.width;
    }    */

    //const cellHeaderTopRight = grid.hitTest(grid.canvas.width-1, 0);
    //const colTopRight = cellHeaderTopRight.gridColumn;
    //const strNameTR = colTopRight.name;

    //const nHColHeader = GridUtils.getColumnHeaderHeight(grid);
    //const cellTopLeft = grid.hitTest(0, nHColHeader);
    const nRowMin = arRowIdxs.length === 0 ? -1 : arRowIdxs[0];//cellTopLeft.gridRow;


    //const nHRow = GridUtils.getRowHeight(grid);
    //const nHCanvas = grid.canvas.height;
    //let nVisibleRowCount = Math.floor((nHCanvas - nHColHeader)/nHRow);

    //let cellBottomLeft = grid.hitTest(0, nVisibleRowCount);//grid.canvas.height-20);
    //if(cellBottomLeft.d === undefined)
      //  throw new Error("Hit test for table row header cell failed.");

    const nRowMax = arRowIdxs.length === 0 ? -2 : arRowIdxs[arRowIdxs.length-1];//nRowMin + nVisibleRowCount;//cellBottomLeft.gridRow;
    //const nRowCount = grid.dataFrame.rowCount;
    //if(nRowMax >= nRowCount)
      //  nRowMax = nRowCount -1;

    arColRowIdxs[0] = arColIdxs.length === 0 ? -1 : arColIdxs[0];//colTopLeft.idx;
    arColRowIdxs[1] = arColIdxs.length === 0 ? -2 : arColIdxs[arColIdxs.length -1];//nColTopRight;
    arColRowIdxs[2] = nRowMin;
    arColRowIdxs[3] = nRowMax;

   // console.log(strNameTL + "   " + strNameTR)
}


GridUtils.initFilterIcons = function(dframe, arFilterClasses)
{
   // ErrorUtils.verifyType(dframe, DG.DataFrame);
    let col = null;
    let filter = null;
    const nColCount = dframe.columns.length;
    for(var nCol = 0; nCol < nColCount; ++nCol)
    {
        col = dframe.columns[nCol];
        try{filter = FiltersUtils.createFilter(col, arFilterClasses);}
        catch(e)
        {
            const cl = arFilterClasses[nCol];
            throw e;
        }

        if(filter === null && col.name.includes("HT Solubility"))
        {
            let sss = 0;
        }

        FiltersUtils.setHasCompatibleFilter(col, filter !== null);
    }
}



GridUtils.initDefaultCellsSizes = function(grid, bCalculateZoom, configGridRenderers)
{
    const nColCount = GridUtils.getColumnCount(grid);

    const nWPrefDefault = grid.dart.DEFAULT_PREFERRED_COLUMN_WIDTH;
    const nHPrefDefault = grid.dart.DEFAULT_PREFERRED_ROW_HEIGHT;

    let nTotalWidthCols = 0;

    let typeSem = null;

    let nWidth = -1;
    let nHTmp = -1;
    let nHeight = -1;

    let renderer = configGridRenderers.getColRowHeaderRenderer();
    if(renderer !== null && render !== undefined) {
        nWidth = renderer.getPreferredCellWidth();
        if (nWidth < grid.dart.MIN_CELL_SIZE)
            nWidth = grid.dart.MIN_CELL_SIZE;
    }

    renderer = configGridRenderers.getRowHeaderRenderer();
    if(renderer !== null && renderer !== undefined) {
        nWidth = Math.max(nWidth, renderer.getPreferredCellWidth());
        if (nWidth < grid.dart.MIN_CELL_SIZE)
            nWidth = grid.dart.MIN_CELL_SIZE;

        nHTmp = renderer.getPreferredCellHeight();
        nHeight = nHTmp;
    }

    let colGrid = grid.columns.byIndex(0);
    GridUtils.setColumnWidth(colGrid, nWidth);

    for(var nCol = 1; nCol < nColCount; ++nCol)
    {
        colGrid = grid.columns.byIndex(nCol);

        if(colGrid.column !== null && (typeSem = GridUtils.getSemType(colGrid.column)) instanceof SemType)
        {
            renderer = configGridRenderers.getColHeaderRenderer(typeSem.constructor);
            nWidth = renderer.getPreferredCellWidth();

            renderer = configGridRenderers.getCellRenderer(typeSem.constructor);
            nWidth = Math.max(nWidth, renderer.getPreferredCellWidth());

            nHTmp = renderer.getPreferredCellHeight();

            //renderer = configGridRenderers.getRowHeaderRenderer();
            //nHTmp = Math.max(nHTmp, renderer.getPreferredCellHeight());
        }
        /*
        else
        {
            renderer = configGridRenderers.getRowHeaderRenderer();

            nWidth =  Math.max(nWidth, renderer.getPreferredCellWidth());
            nHTmp = renderer.getPreferredCellHeight();
        } */
        //nWidth = this.MIN_COLUMN_WIDTH + Math.floor((this.PREFERRED_COLUMN_WIDTH - this.MIN_COLUMN_WIDTH)*fK);

        if (nWidth < grid.dart.MIN_CELL_SIZE)
            nWidth = grid.dart.MIN_CELL_SIZE;

        nTotalWidthCols += nWidth;
        GridUtils.setColumnWidth(colGrid, nWidth);
        GridUtils.setPreferredColumnWidth(colGrid, nWidth);
        GridUtils.setMinimumColumnWidth(colGrid, 0);

        if (nHTmp > nHeight)
            nHeight = nHTmp;
    }


    GridUtils.setPreferredRowHeight(grid, nHeight);
    GridUtils.setMinimumRowHeight(grid, 0.0);

    if(bCalculateZoom)
    {
        let nWVP = grid.canvas.width;
        let nHVP = grid.canvas.height;
        let nSizeMinVP = Math.min(nWVP, nHVP);

        let fK = nSizeMinVP/nTotalWidthCols;
        if(fK < grid.dart.MIN_ZOOM_VALUE)
            fK = grid.dart.MIN_ZOOM_VALUE;

        if(fK > grid.dart.MAX_ZOOM_VALUE)
            fK = grid.dart.MAX_ZOOM_VALUE;

        grid.dart.m_fZoomValue = NaN;//to cheat equals return
        GridUtils.setZoomValue(grid, fK);
        return;
    }
    else grid.dart.m_fZoomValue = 1.0;

    GridUtils.setRowHeight(grid, nHeight);
    /*
    grid.setOptions({
        rowHeight: nHeight
    });*/

    //Repaint
    //colGrid = grid.columns.byIndex(0);
    //my changes colGrid.width = 150;
}


GridUtils.setDefaultColumnSize = function(grid, configGridRenderers) {
    GridUtils.initDefaultCellsSizes(grid, false, configGridRenderers);
}

GridUtils.autoFitColumns = function(grid, configGridRenderers) {
    GridUtils.initDefaultCellsSizes(grid, true, configGridRenderers);
}




/*
GridUtils.findLastSort = function(grid) {
  if(grid === undefined)
      throw new Error("Grid cannot be undefined.");

  let arCols = grid.columns.d.list;

  let column = null;
  let nColCount = arCols.length;
  for(var nCol=0; nCol<nColCount; ++nCol)
  {
      column = grid.columns.byIndex(nCol);
      if(column.m_nSortOption !== undefined && column.m_nSortOption !== null)
      {
          return {grid_column:nCol, sort_option: column.m_nSortOption, sort_ascend: column.m_bSortAscend};
      }
 }

  return null;
}   */

export class SortOptionCtx {
    constructor(nIdxColGrid, nSortOption, bSortAscend) {
        this.m_nIdxColGrid = nIdxColGrid;
        this.m_nSortOption = nSortOption;
        this.m_bSortAscend = bSortAscend;
    }
}





