import { Tutorial } from './tutorial';


/** A collection of tutorials */
export class Track {
  tutorials: Tutorial[];
  name: string;
  helpUrl: string;
  completed: number = 0;

  constructor(name: string, tutorials: Tutorial[], helpUrl: string) {
    this.name = name;
    this.tutorials = tutorials;
    this.helpUrl = helpUrl;
    tutorials.forEach((t) => t.track = this);
  }
}
