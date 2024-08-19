/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {filter} from 'rxjs/operators';
import * as DG from 'datagrok-api/dg';
import {tryParseJson} from "@datagrok-libraries/utils/src/string-utils";
import {FileInfo} from "datagrok-api/dg";

// Power Search: community-curated, template-based, widget-driven search engine

interface Template {
  template: string;
  url: string;
  regexp?: RegExp;   // cached regexp for the template
}

interface Card {
  id: string;
  name: string;
  widget: string;
  templates: Template[];
}

export async function initTemplates(): Promise<void> {

  return; //TODO: bring it back once we can cache app folders. We don't need to load it on each start.

  //let templatesPath = (await _package.getProperties()).get('searchTemplatePaths');
  let templatesPath = 'System:AppData/PowerPack/search-templates';

  async function loadTemplates(): Promise<void> {
    for (let path in templatesPath.split(';')) {
      //if (await grok.dapi.files.exists(path)) {
      if (true) {
        let files = await grok.dapi.files.list(templatesPath, false, null);

        for (let file of files) {
          let s = await file.readAsString();
          let collection: any = tryParseJson(s) ?? [];
          for (let template of collection.templates)
            templates.push(template);
        }
      }
    }
  }

  grok.events.onFileEdited
    // @ts-ignore
    .pipe(filter((f: FileInfo) => f.path.startsWith('PowerPack/search-templates')))
    .subscribe((_) => {
      templates.length = 0;
      loadTemplates().then((_) => grok.shell.info('Search patterns reloaded.'));
    });

  await loadTemplates();
  //console.log(templates);
}

/// Community-curated template collection
export function templatesSearch(s: string, host: HTMLDivElement): void {
  for (let p of templates)
    for (let t of p.templates) {
      let x = <any>t;
      if (x.regexp == null)
        x.regexp = new RegExp(t.template, 'i');
      let matches = x.regexp.exec(s);

      if (matches !== null) {
        console.log(`match! ${p.name}`)

        let widgetProperties: any = {};
        for (let [k, v] of Object.entries(t))
          if (k != 'template' && k != 'regexp') {
            for (let i = 1; i < matches.length; i++)
              v = v.replace('${' + i + '}', matches[i]);

            widgetProperties[k] = v;
          }

        DG.Func.byName(p.widget).apply().then((w: DG.Widget) => {
          w.props.setAll(widgetProperties);
          host.appendChild(w.root);
        });
      }
    }
}

const templates = [
  {
    id: 'foo-test',
    name: 'Foo Test',
    widget: 'kpiWidget',
    templates: [
      {
        template: 'foo ([0-9]+)',
        caption: '${1}'
      }
    ]
  }
]
