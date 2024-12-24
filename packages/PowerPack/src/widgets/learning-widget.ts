/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


interface PlayItem {
  id: string;
  title: string;
  description: string;
  color: string;
}

interface HelpLink {
  title: string;
  url: string;
}

export class LearningWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.box());
    let wikiList = ui.element('ul');
    wikiList.style.columns = '2';
    wikiList.appendChild(ui.div(help.map(renderWiki)));

    let videoList = ui.element('ul');
    videoList.style.columns = '2';
    videoList.style.listStyle = 'none';
    videoList.style.paddingLeft = '0px';
    videoList.appendChild(ui.div(playlists.map(renderPlaylist)));

    let tabs = ui.tabControl();

    tabs.addPane('VIDEO', ()=> ui.panel([videoList]));
    tabs.addPane('WIKI', ()=> ui.panel([wikiList]));
    const daw = DG.Func.find({name: 'demoAppWidget', package: 'Tutorials'})[0];
    const tw = DG.Func.find({name: 'tutorialWidget', package: 'Tutorials'})[0];
    if (daw)
      tabs.addPane('DEMO', () => ui.waitBox(async () => (await daw.apply()).root));
    if (tw)
      tabs.addPane('TUTORIALS', () => ui.waitBox(async () => (await tw.apply()).root));

    this.root.append(tabs.root);

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Learn');
  }
}

function renderPlaylist(p: PlayItem) {
  let url = `https://www.youtube.com/playlist?list=${p.id}`;
  let listItem = ui.element('li');
  listItem.style.breakInside = 'avoid';
  listItem.style.pageBreakInside = 'avoid';

  let img = ui.iconFA('play-circle');
  img.style.backgroundColor = p.color;
  img.classList.remove('fal');
  img.classList.add('far');
  img.classList.add('power-pack-learn-icon');

  let link = ui.link(p.title,url,p.description,'');
  let root = ui.divH([img, link], {style: { alignItems: 'center'}});

  listItem.appendChild(root);
  return listItem;
}

function renderWiki (p: HelpLink) {
  let listItem = ui.element('li');
  listItem.appendChild(ui.link(p.title, p.url));
  return listItem;
}

// get info on channels:
// https://developers.google.com/youtube/v3/docs/playlists/list?apix_params=%7B%22part%22%3A%5B%22snippet%22%5D%2C%22channelId%22%3A%22UCXPHEjOd4gyZ6m6Ji-iOBYg%22%7D

let playlists: PlayItem[] = [
  {
    id: "PLIRnAn2pMh3kvsE5apYXqX0I9bk257_eY",
    title: "Meetings",
    description: "Join us for the regular community meetings",
    color: '#c9b6b6'
  },
  {
    id: "PLIRnAn2pMh3ncmRaE4fjmPbDaDvKklh7j",
    title: "Develop",
    description: "Extend the platform by developing scripts, functions, and packages",
    color: "#cec3ad"
  },
  {
    id: "PLIRnAn2pMh3nHUxed6p-uw7If24nGENDa",
    title: "Cheminformatics",
    description: "Work with chemical structures using the first-class cheminformatics support provided by the platform",
    color: "#a8c0bf"
  },
  {
    id: "PLIRnAn2pMh3nLvDs3NLXkLtsyJeX912GG",
    title: "Visualize",
    description: "Understand your data by using powerful interactive visualizations",
    color: "#a2acc7"
  },
  {
    id: "PLIRnAn2pMh3nToHhFs3eXpf9xXa195lrN",
    title: "Explore",
    description: "Discover, transform, and explore your data, no matter where it comes from",
    color: '#b8b0cd'
  }
];



let help: HelpLink[] = [
  {
    title: 'Cheminformatics',
    url: 'https://datagrok.ai/help/datagrok/solutions/domains/chem/'
  },
  {
    title: 'Bioinformatics',
    url: 'https://datagrok.ai/help/datagrok/solutions/domains/bio/'
  },
  {
    title: 'Overview',
    url: 'https://datagrok.ai/help/datagrok/'
  },
  {
    title: 'Access',
    url: 'https://datagrok.ai/help/access/'
  },
  {
    title: 'Transform',
    url: 'https://datagrok.ai/help/transform/'
  },
  {
    title: 'Visualize',
    url: 'https://datagrok.ai/help/visualize/viewers/'
  },
  {
    title: 'Explore',
    url: 'https://datagrok.ai/help/explore/cluster-data'
  },
  {
    title:'Compute',
    url: 'https://datagrok.ai/help/compute/'
  },
  {
    title: 'Predict',
    url: 'https://datagrok.ai/help/learn/'
  },
  {
    title: 'Discover',
    url: 'https://datagrok.ai/help/explore/data-augmentation/'
  },
  {
    title: 'Develop',
    url: 'https://datagrok.ai/help/develop/'
  },
  {
    title: 'Govern',
    url: 'https://datagrok.ai/help/govern/audit'
  },
  {
    title: 'Acknowledgments',
    url: 'https://datagrok.ai/help/datagrok/resources/acknowledgements'
  }
];
