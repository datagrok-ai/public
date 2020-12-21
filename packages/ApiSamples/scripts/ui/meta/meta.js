// Custom meta classes
// Clicking on the item shows its properties in the property panel
// Right-clicking allows to use context actions

class Fruit {
  constructor(name) {
    this.name = name;
  }
}

// Defines the way Datagrok handles entities of the specified type
class FruitMeta extends DG.JsEntityMeta {
  get type() {
    return 'fruit'
  }

  // Checks whether this meta class is the handler for [x]
  isApplicable(x) {
    return x instanceof Fruit;
  }

  renderCard(x) {
    return ui.bind(x, ui.divText(`Fruit: ${x.name}`));
  }

  // Renders properties in the property panel
  renderProperties(x) {
    return ui.divText(`Properties for ${x.name}`);
  }

  // Renders tooltip
  renderTooltip(x) {
    return ui.divText(`${x.name} is in the air!`);
  }

  init() {
    this.registerParamFunc('Eat', (fruit) => {
      grok.shell.info(`Ate ${fruit.name}`);
    });
  }
}

// Register meta class with the platform
DG.JsEntityMeta.register(new FruitMeta());

grok.shell.newView().append(ui.div([
  ui.renderCard(new Fruit('apple')),
  ui.renderCard(new Fruit('banana')),
]));
