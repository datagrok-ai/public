import * as DG from 'datagrok-api/dg';


/** A base class for tutorials */
export abstract class Tutorial extends DG.Widget {
  abstract get name(): string;
  abstract get description(): string;

  demoTable: string = 'demog';
}


/** A collection of tutorials */
export class Track {
  tutorials: Tutorial[];
  name: string;

  constructor(name: string, ...tutorials: Tutorial[]) {
    this.name = name;
    this.tutorials = tutorials;
  }
}
