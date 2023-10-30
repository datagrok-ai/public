//name: cars
//description: Demonstrates value lookup
//language: javascript
//input: string model { choices: OpenFile("System:AppData/Samples/cars.csv"); propagateChoice: all }
//input: double mpg
//input: int cyl
//input: int disp
//input: int hp
//input: double drat
//input: double wt
//input: double qsec
//output: double mpgPerCylinder
//output: double acc

mpgPerCylinder = mpg / cyl;
acc = hp / wt;