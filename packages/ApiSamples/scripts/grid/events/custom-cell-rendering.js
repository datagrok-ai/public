let view = grok.shell.addTableView(grok.data.demo.demog());

// custom Canvas-based rendering

view.grid.onCellRender.subscribe(function (args) {
  args.g.beginPath();
  args.g.arc(args.bounds.x + args.bounds.width / 2, args.bounds.y + args.bounds.height / 2, 10, 0, Math.PI * 2, true);
  args.g.closePath();
  args.g.fillStyle = (args.cell.isColHeader ? 'red' : 'blue');
  args.g.fill();
  args.preventDefault();
});
