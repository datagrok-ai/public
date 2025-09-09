import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import {div} from "datagrok-api/ui";
import {IDartApi} from "datagrok-api/src/api/grok_api.g";

const api: IDartApi = <any>window;

export class RecentProjectsWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.panel());
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent projects');
    this.getProjects();
  }

  async getProjects(){
    const projectList = ui.wait(async () => {
      const listDiv = ui.divV([]);
      const projects = await grok.dapi.projects.recent.list({pageSize: 5});
      if (projects.length > 0)
        listDiv.appendChild(ui.divV(projects.map((p) => getProjectCard(p) ), {style:{overflowX:'scroll'}}));
      else {
        this.caption = 'Demo projects';
        grok.dapi.projects
          .filter('#demo')
          .list({pageSize: 5})
          .then((projects) => listDiv.appendChild(ui.divV(projects.map((p) => getDemoProjectCard(p)), {style: {overflowX: 'scroll'}})));
      }
      if (projects.length < 20) {
        const dropZone = div([
          ui.link('Open a file', () => grok.shell.openFileOpenDialog()),
          ui.span([', or drop it here'])
        ], {classes: 'pp-drop-zone'});
        listDiv.appendChild(dropZone);
      }
      return listDiv;
    });
    this.root.appendChild(projectList);
  }
}

function getProjectCard(p: DG.Project): HTMLElement {
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, (_: any) => {}, undefined, null))
    ]);
  ui.bind(p, card, {contextMenu: true});
  return card;
}

function getDemoProjectCard(p: DG.Project): HTMLElement {
  let card = ui.cards.summary(
    ui.image(p.pictureUrl, 70, 50, {target: () => {}}),
    [
      ui.h2(ui.link(p.friendlyName, p, undefined, null)),
      ui.div([p.createdOn], {style:{color:'var(--grey-4)'}})
    ]);
  ui.bind(p, card);
  return card;
}
