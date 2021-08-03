//Example of various rendering items

let view = grok.shell.newView('Rendering items');

let table = grok.data.testData('demog');
let grid1 = DG.Viewer.grid(table);
let grid2 = DG.Viewer.grid(table);

// Grid render
grid1.onCellRender.subscribe(function (args) {
  args.g.beginPath();
  args.g.arc(args.bounds.x + args.bounds.width / 2, args.bounds.y + args.bounds.height / 2, 10, 0, Math.PI * 2, true);
  args.g.closePath();
  args.g.fillStyle = (args.cell.isColHeader ? 'red' : 'blue');
  args.g.fill();
  args.preventDefault();
});

grid2.onCellRender.subscribe(function (args) {
  args.g.fillStyle = (args.cell.isColHeader ? 'red' : (args.cell.isRowHeader ? 'green' : 'blue'));
  args.g.fillText(`${args.cell.tableRowIndex}`, args.bounds.x + args.bounds.width / 2, args.bounds.y + args.bounds.height / 2);
  args.preventDefault();
});

// Render card
let rendedCard = ui.div();
grok.dapi.scripts
  .list({pageSize: 1})
  .then((scripts) => rendedCard.append(ui.div(scripts.map((p) => ui.renderCard(p)))));

// Render list
let rendedList = ui.div();
grok.dapi.users
  .list({pageSize: 5})
  .then((users) => rendedList.append(ui.list(users)))

// Render table  
let rendedTable = ui.div();
grok.dapi.users.current().then((user) => rendedTable.append(ui.tableFromMap({
  'First name': user.firstName,
  'Last name': user.lastName,
  'Login': user.login,
  'Picture': user.picture,
  'Project': user.project,
})));

// Layout
let container = ui.block([
  ui.block50([
    ui.h3('Canvas-based rendering'),
    grid1.root
  ]),
  ui.block50([
    ui.h3('Canvas-based rendering'),
    grid2.root
  ]),
  ui.block50([
    ui.h3('Molecule render'),
    grok.chem.svgMol('c1(ccc2N=C(C)N(C(=O)c2c1)c3ccc(OC)cc3)NC(=S)Nc4ccccc4')
  ]),
  ui.block50([
    ui.h3('Render card'),
    rendedCard
  ]),
  ui.block50([
    ui.h3('Render list'),
    rendedList
  ]),
  ui.block50([
    ui.h3('Render table'),
    rendedTable,
  ]),
]);
view.append(container);