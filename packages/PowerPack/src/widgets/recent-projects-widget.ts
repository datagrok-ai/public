import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";

export class RecentProjectsWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.div());

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Recent projects');

    let projectList = ui.divV([]);
    grok.dapi.projects
      .filter('#demo')
      .list({pageSize: 5})
      .then((projects) => projectList.appendChild(ui.divV(projects.map((p) => projectData(p) ), {style:{overflowX:'scroll'}})));

    this.root.appendChild(projectList);
  }
}

function projectData(p: DG.Project){
  let root = ui.divV([]);
  let ppath = p.path;
  //let pictureUrl = 'https://dev.datagrok.ai/api/entities/picture/'+p.d.PictureMixin_pictureId+'?w=200';
  let picture = ui.element('img') as HTMLImageElement;
  let name = p.friendlyName;
  let description = p.d.AuthorMixin_createdOn;
  picture.src = p.pictureUrl;
  picture.width = 70;
  picture.height  = 50;
  picture.style.marginRight = '10px';
  picture.style.border = '1px solid var(--grey-1)';
  picture.style.borderRadius = '1px';

  root.style.cursor='pointer';
  root.style.flex='none';
  root.style.margin="0";
  root.style.padding="5px";
  root.setAttribute('onmouseover','this.style.background="#f5f5f5"');
  root.setAttribute('onmouseout','this.style.background="none"');
  root.appendChild(
    ui.divH([
      picture,
      ui.divV([
        ui.div(name, {style:{color:'var(--blue-1)',marginBottom:'5px'}}),
        ui.div(description, {style:{color:'var(--grey-4)'}})
      ], {style:{justifyContent:'center'}})
    ])
  );

  ui.bind(p, root);
  // ui.tooltip.bind(root, ui.divV([ui.span(['Name: ',name]),ui.span(['Description: ',p.description]),ui.span(['User: ',p.d.AuthorMixin_author]),ui.span(['Created: ',p.d.AuthorMixin_createdOn])])),
  //   root.onclick = (e) =>window.location.href='https://dev.datagrok.ai'+ppath;

  return root;
}