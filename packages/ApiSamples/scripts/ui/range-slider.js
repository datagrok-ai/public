// Creates a range slider, subscribes to its stream, which fires new events when values change

function info(s) { grok.shell.info(s); }

let v = grok.shell.newView('demo: range slider');

let rangeSlider = ui.rangeSlider(0, 10, 2, 5);
rangeSlider.onValuesChanged.subscribe((_) => info('rangeSlider-values-changed'));

v.append(rangeSlider);