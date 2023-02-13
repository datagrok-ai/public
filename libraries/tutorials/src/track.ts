import {Tutorial} from './tutorial';


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

  async updateStatus(): Promise<{[id: string]: boolean}> {
    const statusMap: {[id: string]: boolean} = {};
    let completed = 0;
    for (let i = 0; i < this.tutorials.length; i++) {
      const t = this.tutorials[i];
      await t.updateStatus();
      statusMap[i] = t.status;
      if (t.status)
        completed++;
    }
    this.completed = completed;
    return statusMap;
  }
}
