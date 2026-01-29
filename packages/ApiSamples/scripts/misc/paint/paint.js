// Helper methods for rendering on canvas

const canvas = ui.canvas(200, 200);

const g = canvas.getContext('2d');
DG.Paint.horzAxis(g, 0, 100, 30, 170, 160, 30);
DG.Paint.vertAxis(g, 0, 100, 0, 10, 30, 160);

for (let i = 0; i < 100; i++) {
  DG.Paint.marker(g, DG.MARKER_TYPE.CIRCLE,
    Math.random() * 160 + 30, Math.random() * 160 + 10,
    DG.Color.filteredRows, 5);
}

canvas;