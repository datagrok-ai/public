// This examples shows how to register a function that becomes a first-class
// citizen in the platform (i.e., it can be used from console, gets registered
// in help, there could be an optional audit trail associated with the invocations, etc)
//
// The code below registers two functions, "jsConcat" and "jsWidget". To test
// jsConcat, enter "jsConcat(42, 33)" in the console(https://datagrok.ai/help/overview/navigation#console).
// To test jsWidget, create a new Dashboard, and click on "Widget" under "Widgets".
//
// See also https://datagrok.ai/help/overview/functions/function,
// https://datagrok.ai/help/develop/js-api#registering-functions

grok.functions.register({
  signature: 'String jsConcat(int foo, int bar)',
  run: (foo, bar) => `${foo}_${bar}`
});

grok.functions.register({
  signature: 'double norm2(double x, double y)',
  run: (x, y) => Math.sqrt(x / y)
});

grok.functions.register({
  signature: 'String/Molecule testMol()',
  run: () => `C(=O)(O)c1ccccc1N`
});

grok.functions.register({
  signature: 'Widget jsWidget()',
  tags: 'Widgets',
  run: function () {
    let e = document.createElement('DIV');

    function update() {
      let date = new Date();
      e.innerText = date.toTimeString();
    }

    window.setTimeout(update, 1000);

    return new DG.Widget(e);
  }
});
