import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class RecentProjectsWidget extends DG.Widget {
  caption: string;
  order: string;

  constructor() {
    super(ui.panel());
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent projects');
    this.order = super.addProperty('order', DG.TYPE.STRING, '1');
    this.getProjects();
  }

  async getProjects(){
    let count = (await grok.dapi.projects.recent.list()).length;
    let projectList = ui.divV([]); 

    if (count != 0){
      grok.dapi.projects.recent
      .list({pageSize: 5})
      .then((projects) => projectList.appendChild(ui.divV(projects.map((p) => projectData(p) ), {style:{overflowX:'scroll'}})));
    } else {
      this.caption = 'Demo projects';
      grok.dapi.projects
      .filter('#demo')
      .list({pageSize: 5})
      .then((projects) => projectList.appendChild(ui.divV(projects.map((p) => demoProjectData(p) ), {style:{overflowX:'scroll'}})));
    }

    this.root.appendChild(projectList);
  }
}

function projectData(p: DG.Project){
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, (_: any) => {}, undefined, null))
    ]);
  ui.bind(p, card, {contextMenu: true});
  return card;
}

function demoProjectData(p: DG.Project){
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, p, undefined, null)),
      ui.div([p.createdOn], {style:{color:'var(--grey-4)'}})
    ]);
  ui.bind(p, card);
  return card;
}
