/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class LearningWidget extends DG.Widget {
  caption: string;
  order: string;

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

    let tabs = ui.tabControl({
       'VIDEO':ui.panel([
        videoList,
       ]),
       'WIKI':ui.panel([
        wikiList
       ])
     });
    
    grok.functions
    .call("Tutorials:tutorialWidget")
    .then((res) => {
      tabs.addPane('TUTORIALS',()=>res._root)
    })
    .catch((error) => {
      console.log('TutorialWidget not found. \n'+error);
    });
    
    this.root.append(tabs.root);

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Learn');
    this.order = super.addProperty('order', DG.TYPE.STRING, '5');
  }
}

function renderPlaylist(p: any) {
  let url = `https://www.youtube.com/playlist?list=${p.id}`;
  let listItem = ui.element('li');
  listItem.style.breakInside = 'avoid';
  listItem.style.pageBreakInside = 'avoid';
  listItem.style.minWidth = '155px';

  let root = ui.divH([]);
  root.style.alignItems = 'center';

  let img = ui.image(p.url, 35,35,{target:url});
  img.style.marginRight = '10px';

  let link = ui.link(p.title,url,p.description,'');

  root.appendChild(img);
  root.appendChild(link);
  listItem.appendChild(root);
  return listItem
  // return ui.cards.summary(
  //   ui.image(p.url, 120, 90, { target: url }),
  //   [
  //     ui.h2(ui.link(p.title, url)),
  //     ui.divText(p.description)
  // ]);
}

function renderWiki (p: any){
  let listItem = ui.element('li');
  listItem.appendChild(ui.link(p.title, p.url));
  return listItem;
}

// get info on channels:
// https://developers.google.com/youtube/v3/docs/playlists/list?apix_params=%7B%22part%22%3A%5B%22snippet%22%5D%2C%22channelId%22%3A%22UCXPHEjOd4gyZ6m6Ji-iOBYg%22%7D

let playlists = [
  {
    "id": "PLIRnAn2pMh3kvsE5apYXqX0I9bk257_eY",
    "title": "Meetings",
    "description": "Join us for the regular community meetings",
    "url": "https://raw.githubusercontent.com/datagrok-ai/public/master/help/uploads/pictures/meetings-th.png"
  },
  {
    "id": "PLIRnAn2pMh3ncmRaE4fjmPbDaDvKklh7j",
    "title": "Develop",
    "description": "Extend the platform by developing scripts, functions, and packages",
    "url": "https://raw.githubusercontent.com/datagrok-ai/public/master/help/uploads/pictures/develop-th.png"
  },
  {
    "id": "PLIRnAn2pMh3nHUxed6p-uw7If24nGENDa",
    "title": "Cheminformatics",
    "description": "Work with chemical structures using the first-class cheminformatics support provided by the platform",
    "url": "https://raw.githubusercontent.com/datagrok-ai/public/master/help/uploads/pictures/chem-th.png"
  },
  {
    "id": "PLIRnAn2pMh3nLvDs3NLXkLtsyJeX912GG",
    "title": "Visualize",
    "description": "Understand your data by using powerful interactive visualizations",
    "url": "https://raw.githubusercontent.com/datagrok-ai/public/master/help/uploads/pictures/visualize-th.png"
  },
  {
    "id": "PLIRnAn2pMh3nToHhFs3eXpf9xXa195lrN",
    "title": "Explore",
    "description": "Discover, transform, and explore your data, no matter where it comes from",
    "url": "https://raw.githubusercontent.com/datagrok-ai/public/master/help/uploads/pictures/explore-th.png"
  }
];



let help = [
  {
    'title': 'Access',
    'url': 'https://datagrok.ai/help/access/db-exploration'
  },
  {
    'title':'Collaborate',
    'url': 'https://datagrok.ai/help/collaborate/chat'
  },
  {
    'title': 'Develop',
    'url': 'https://datagrok.ai/help/develop/getting-started'
  },
  {
    'title': 'Discover',
    'url': 'https://datagrok.ai/help/discover/data-augmentation'
  },
  {
    'title': 'Domains',
    'url': 'https://datagrok.ai/help/domains/chem/chemically-aware-viewers'
  },
  {
    'title': 'Explore',
    'url': 'https://datagrok.ai/help/explore/cluster-data'
  },
  {
    'title': 'Govern',
    'url': 'https://datagrok.ai/help/govern/audit'
  },
  {
    'title': 'Learn',
    'url': 'https://datagrok.ai/help/learn/data-science'
  },
  {
    'title': 'Overview',
    'url': 'https://datagrok.ai/help/datagrok/create-project'
  },
  {
    'title': 'Stories',
    'url': 'https://datagrok.ai/help/solutions/life-sciences'
  },
  {
    'title': 'Transform',
    'url': 'https://datagrok.ai/help/transform/add-new-column'
  },
  {
    'title': 'Visualize',
    'url': 'https://datagrok.ai/help/visualize/add-custom-form'
  },
  {
    'title': 'Acknowledgments',
    'url': 'https://datagrok.ai/help/acknowledgements'
  }
];

