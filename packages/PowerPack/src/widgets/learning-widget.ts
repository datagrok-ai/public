/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class LearningWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.panel());

    this.root.appendChild(ui.divV(playlists.map(renderPlaylist)));

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Learn');
  }
}

function renderPlaylist(p: any) {
  let url = `https://www.youtube.com/playlist?list=${p.id}`;
  return ui.cards.summary(
    ui.image(p.url, 120, 80, { target: url }),
    [
      ui.h2(ui.link(p.title, url)),
      ui.divText(p.description)
  ]);
}

// get info on channels:
// https://developers.google.com/youtube/v3/docs/playlists/list?apix_params=%7B%22part%22%3A%5B%22snippet%22%5D%2C%22channelId%22%3A%22UCXPHEjOd4gyZ6m6Ji-iOBYg%22%7D

let playlists = [
  {
    "id": "PLIRnAn2pMh3kvsE5apYXqX0I9bk257_eY",
    "title": "Meetings",
    "description": "Join us for the regular community meetings",
    "url": "https://i.ytimg.com/vi/JaJgxtHAb98/default.jpg"
  },
  {
    "id": "PLIRnAn2pMh3ncmRaE4fjmPbDaDvKklh7j",
    "title": "Develop",
    "description": "Extend the platform by developing scripts, functions, and packages",
    "url": "https://i.ytimg.com/vi/dKrCk38A1m8/default.jpg"
  },
  {
    "id": "PLIRnAn2pMh3nHUxed6p-uw7If24nGENDa",
    "title": "Cheminformatics",
    "description": "Work with chemical structures using the first-class cheminformatics support provided by the platform",
    "url": "https://i.ytimg.com/vi/wCdzD64plEo/default.jpg"
  },
  {
    "id": "PLIRnAn2pMh3nLvDs3NLXkLtsyJeX912GG",
    "title": "Visualize",
    "description": "Understand your data by using powerful interactive visualizations",
    "url": "https://i.ytimg.com/vi/jHRpOnhBAz4/default.jpg"
  },
  {
    "id": "PLIRnAn2pMh3nToHhFs3eXpf9xXa195lrN",
    "title": "Explore",
    "description": "Discover, transform, and explore your data, no matter where it comes from",
    "url": "https://i.ytimg.com/vi/67LzPsdNrEc/default.jpg"
  }
];

