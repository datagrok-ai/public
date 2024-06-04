// Solver's callback base

export class Callback {
  private count = 1;
  constructor() {};

  public onIterationStart(): void {
    console.log(`Call No. ${this.count}`);
    this.count++;
  }
};
