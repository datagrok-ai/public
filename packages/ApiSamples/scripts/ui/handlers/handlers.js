// Custom object handlers
// Clicking on the item shows its properties in the context panel
// Right-clicking allows to use context actions

class Fruit {
  constructor(name, color) {
    this.name = name;
    this.color = color;
  }
}

class FruitCanvasRenderer extends DG.CanvasRenderer {
  render(g, x, y, w, h, fruit, context) {
    g.fillText(fruit.name, x + 10, y + 10);
  }
}

class FruitGridCellRenderer extends DG.GridCellRenderer {
  get cellType() {return 'fruit';}
  render(g, x, y, w, h, cell, style) {
    g.fillStyle = cell.cell.value.color;
    g.fillText(cell.cell.value.name, x + 10, y + 10);
  }
}

// Defines the way Datagrok handles entities of the specified type
class FruitHandler extends DG.ObjectHandler {
  get type() {return 'fruit';}

  // Checks whether this is the handler for [x]
  isApplicable(x) {return x instanceof Fruit;}

  getCanvasRenderer(x) {return new FruitCanvasRenderer();}
  getGridCellRenderer(x) {return new FruitGridCellRenderer();}

  renderIcon(x) {return ui.iconFA('apple-alt');}
  renderMarkup(x) {let m = ui.span([this.renderIcon(), ui.label(x.name)]); $(m).css('color', x.color); return m;}
  renderProperties(x) {return ui.divText(`Properties for ${x.name}`);}
  renderTooltip(x) {return ui.divText(`${x.name} is in the air!`);}
  renderCard(x, context) {
    return ui.bind(x, ui.divV([
      this.renderMarkup(x),
      ui.divText(`Context: ${context}`)
    ]), 'd4-gallery-item');
  }

  init() {
    this.registerParamFunc('Eat', (fruit) => {
      grok.shell.info(`Ate ${fruit.name}`);
    });
  }
}

// Register handler with the platform
DG.ObjectHandler.register(new FruitHandler());

let apple = new Fruit('apple', 'red');
let orange = new Fruit('orange', 'orange');
let v = grok.shell.newView();

// automatic rendering
v.append(ui.div([
  ui.h1('Automatic rendering'),
  ui.renderCard(apple),
  ui.renderCard(orange),
]));

// manual rendering
let handler = DG.ObjectHandler.forEntity(orange);
v.append(ui.div([
  ui.h1('Manual rendering'),
  ui.tableFromMap({
    'icon': handler.renderIcon(orange),
    'markup': handler.renderMarkup(orange, 'myContext'),
    'card': handler.renderCard(orange),
    'tooltip': handler.renderTooltip(orange),
    'properties': handler.renderProperties(orange),
  })
]));

v.append(ui.div([
  ui.h1('Inline rendering'),
  ui.inlineText(['I like to eat ', orange, ' in the morning'])
]));

let canvas = ui.canvas(200, 100);
handler.getCanvasRenderer().render(canvas.getContext('2d'), 0, 0, 200, 100, apple, null);
v.append(ui.div([
  ui.h1('Canvas rendering'),
  canvas
]));

// A column gets associated with the handler via semType
let fruitColumn = DG.Column.fromType(DG.TYPE.OBJECT, 'fruits', 2);
fruitColumn.init((i) => [apple, orange][i]);
fruitColumn.semType = 'fruit';
let table = DG.DataFrame.fromColumns([fruitColumn]);
v.append(ui.h1('Grid'));
v.append(DG.Viewer.grid(table).root);