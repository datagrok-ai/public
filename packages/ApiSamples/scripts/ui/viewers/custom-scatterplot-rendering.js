// Renders additional lines that connect points

function renderLines (sp) {
    let ctx = sp.getInfo()['canvas'].getContext('2d');
    ctx.beginPath();
    ctx.strokeStyle = 'black';

    for (let i = 0; i < sp.dataFrame.rowCount; i++) {
        let point = sp.coordsToScreen(sp.dataFrame.get('#0', i), sp.dataFrame.get('#1', i));
        ctx.lineTo(point.x, point.y);
    }

    ctx.stroke();
}

let sp = DG.Viewer.scatterPlot(grok.data.demo.randomWalk(20, 2), {
    xColumnName: '#0',
    yColumnName: '#1'
});
grok.shell.newView('View', [sp]);

sp.onEvent('d4-after-draw-scene').subscribe(_ => renderLines(sp));