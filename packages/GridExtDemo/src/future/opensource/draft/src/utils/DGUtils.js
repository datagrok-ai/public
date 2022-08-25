
export class DGUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}

DGUtils.isCurrentView = function(view)
{
    const b = grok.shell.v.dart === view.dart;
    return b;
}

DGUtils.containsViewer = function(view, viewer)
{
    const lstViewers = view.viewers;
    for(var viewerTmp of lstViewers)
    {
        if(viewerTmp.dart === viewer.dart)
            return true;
    }

    return false;
}


DGUtils.getViewerTypeInstanceCount = function(view, strType)
{
    let nCount = 0;
    const lstViewers = view.viewers;
    for(var viewerTmp of lstViewers)
    {
        if(viewerTmp.type === strType)
            ++nCount;
    }

    return nCount;
}


DGUtils.findViewersbyType = function(view, strType)
{
   // ErrorUtils.verifyType(view, DG.TableView)

    let ar = [];
    const lstViewers = view.viewers;
    for(var viewer of lstViewers)
    {
        if(viewer.type === strType)
            ar.push(viewer);
    }

    return ar;
}


DGUtils.iterateDockNodes = function(node, set, nCount)
{
    if(nCount === undefined)
        nCount = 0;
    else ++nCount;

    const c = node.container;
    const e = c.containerElement;
    const lstViewerRoots = e.querySelectorAll('.d4-viewer');
    const strName = e.nodeName;
    console.log(nCount + " " + strName + " " + e.className);
    for(let nE=0; nE<lstViewerRoots.length; ++nE) {
     set.add(lstViewerRoots.item(nE));
    }

    const itChildren = node.children;
    const arChildren = Array.from(itChildren);
    let nodeChild = null;
    for(let n=0; n<arChildren.length; ++n) {
        nodeChild = arChildren[n];
        nCount = DGUtils.iterateDockNodes(nodeChild, set, nCount);
    }


    return nCount;
}
