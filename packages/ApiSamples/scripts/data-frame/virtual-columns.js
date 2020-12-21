class Car {
  constructor(make, model, price) {
    this.make = make;
    this.model = model;
    this.price = price;
    console.log('car constructed');
  }

  toString() {
    return `${this.make} ${this.model}: ${this.price}`;
  }
}

let table = DG.DataFrame.fromCsv(
  `make, model,    cylinders, volume, price
Honda, Civic,    4,         1.4,    15000
Honda, Accord,   6,         1.8,    20000
BMW,   328i,     4,         1.7,    60000        
BMW,   535i,     6,         1.5,    35000
Tesla, Roadster, ,          1.6,    100000
Tesla, Model S,  ,          1.6,    120000`);

table.columns.addNewVirtual('car', (i) => new Car(table.row(i).make, table.row(i).model, table.row(i).price));

grok.shell.add(table);