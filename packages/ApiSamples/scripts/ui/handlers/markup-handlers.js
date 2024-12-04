//name: markup-handlers
//tags: demo
//language: javascript
// Custom markup handlers
// Allows logging and audit of events linked to in-app objects

// Defines a sample class
class Fruit {
  constructor(name, color) {
    this.name = name;
    this.color = color;
  }
}

// Defines the way Datagrok handles entities of the specified type
// Shortened version of a sample handlers.js
class FruitHandler extends DG.ObjectHandler {
  get type() { return 'fruit' }

  // Checks whether this is the handler for [x]
  isApplicable(x) { return x instanceof Fruit; }

  renderIcon(x) { return ui.iconFA(`apple-alt`); }
  renderMarkup(x) { let m = ui.span([this.renderIcon(x), ui.label(x.name)]); $(m).css('color', x.color); return m; }

  // Is used by platform to detect markup description of Fruit
  get markupRegexp() {
    return "fruit\\.([a-zA-Z]+)\\.([a-zA-Z]+)";
  }

  // Converts fruit to markup description
  toMarkup(x) {
     return `fruit.${x.name}.${x.color}`;
  }

  // Deserialize object from parsed markup
  // Accepts list of strings, that comes from regexp group matches
  // In our case, matches = ["fruit.apple.red", "apple", "red"]
  fromMarkup(matches) {
    var name = matches[1];
    var color = matches[2];
    return new Fruit(name, color);
  }
}

// Register handler with the platform
var handler = new FruitHandler();
DG.ObjectHandler.register(handler);

// Create some fruits
let apple = new Fruit('apple', 'red');
let orange = new Fruit('orange', 'orange');

// Log usage of the fruit.
// Highlight custom markup objects with "#"
grok.log.usage('@user liked #fruit', {'user': DG.User.current(), '#fruit': apple});
grok.log.usage('@user disliked #fruit', {'user': DG.User.current(), '#fruit': orange});

let logObject = await grok.dapi.log.filter('description = "@user liked #fruit"').first();
grok.shell.info(ui.render(logObject));