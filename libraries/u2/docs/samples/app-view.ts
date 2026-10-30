import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AsyncValue, CatalogTable, TextInput, computed, iconButton} from '@datagrok-libraries/u2';
import {ViewHost} from '@datagrok-libraries/u2/dg';

//name: Projects
//tags: app
//meta.browsePath: Demo
//output: view result
export function projectsApp(): DG.ViewBase {
  return ViewHost.create({name: 'Projects', content: () => {
    const search = new TextInput({search: true, inline: true, placeholder: 'Search projects'});
    const projects = new AsyncValue(() => grok.dapi.projects.filter(search.value).list({pageSize: 200}),
      {deps: () => search.value});
    const refresh = iconButton('sync', () => projects.refresh(), {tooltip: 'Refresh'});
    const table = new CatalogTable<DG.Project>({
      name: 'projects', items: projects, keyOf: (p) => p.id,
      columns: [
        {name: 'Name', value: (p) => p.friendlyName, pinned: true, link: 'Open the project'},
        {name: 'Author', value: (p) => p.author?.friendlyName},
        {name: 'Created', value: (p) => p.createdOn.toDate()},
        {name: 'Description', value: (p) => p.description, maxWidth: '360px'},
      ],
      groupings: [{name: 'Author', key: (p) => p.author?.friendlyName}],
      onActivate: (p) => p.open(),
      persist: 'my-app.projects',
    });
    return {root: table, ribbon: [[search], [refresh]], status: computed(() => `${projects.value?.length ?? '…'} projects`),
      path: '/apps/MyApp/Projects'};
  }}).view;
}
