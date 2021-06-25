import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class RecentProjectsWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.panel());

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent projects');

    let projectList = ui.divV([]); 

    let projectsCount = grok.dapi.projects.recent.include.length;

    if (projectsCount != 0){
      grok.dapi.projects.recent
      .list({pageSize: 5})
      .then((projects) => projectList.appendChild(ui.divV(projects.map((p) => projectData(p) ), {style:{overflowX:'scroll'}})));
    } else {
      this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Demo projects');
      grok.dapi.projects
      .filter('#demo')
      .list({pageSize: 5})
      .then((projects) => projectList.appendChild(ui.divV(projects.map((p) => demoProjectData(p) ), {style:{overflowX:'scroll'}})));
    }

    this.root.appendChild(projectList);
  }
}

function projectData(h: DG.HistoryEntry){
  let p = <DG.Project>h.object;
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, p)),
      ui.div([h.time], {style:{color:'var(--grey-4)'}})
    ]);
  ui.bind(p, card);
  return card;
}

function demoProjectData(p: DG.Project){
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, p)),
      ui.div([p.createdOn], {style:{color:'var(--grey-4)'}})
    ]);
  ui.bind(p, card);
  return card;
}