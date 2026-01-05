const properties = [
  DG.Property.jsInt('age'),
  DG.Property.jsString('sex'),
  DG.Property.jsBool('control')
];

const items = [
  {age: 28, sex: 'M', control: true},
  {age: 35, sex: 'M', control: false},
  {age: 56, sex: 'F', control: true},
  {age: 30, sex: 'F', control: false}
];

DG.Grid.fromProperties(items, properties).root;