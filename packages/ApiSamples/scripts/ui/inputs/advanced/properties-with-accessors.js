// Properties with custom accessors: get/set win over the default field read/write

const object = {celsius: 20};

const properties = [
  DG.Property.js('celsius', DG.TYPE.FLOAT, {units: '°C'}),
  DG.Property.js('fahrenheit', DG.TYPE.FLOAT, {
    units: '°F',
    get: (o) => o.celsius * 9 / 5 + 32,
    set: (o, v) => o.celsius = (v - 32) * 5 / 9,
  }),
];

ui.input.form(object, properties);
