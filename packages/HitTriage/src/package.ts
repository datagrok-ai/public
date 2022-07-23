/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Template} from "./template";
import {TemplateHandler} from "./template_handler";
import {hitTriageView} from "./hit-triage-view";

export let _package = new DG.Package();
export let templates: Array<Template>;

//tags: init
export function init() {
  DG.ObjectHandler.register(new TemplateHandler());
  grok.functions.onAfterRunAction.subscribe((call) => {
    let template = call.context.getVariable('template');
    if (template == null)
      return;
    console.log(call.func.outputs);
    call.func.outputs.forEach((param)  => {
      let output = param.get(call);
      if (output instanceof DG.DataFrame)
        output.setTag('HitTriage', template.name);
    });
  });
}

//tags: app
//name: Hit Triage
export async function hitTriageApp() {
  grok.shell.addView(hitTriageView());
  // let v = grok.shell.newView('Hit Triage Templates');
  // let templateFiles = await _package.files.list('templates/', false);
  // templates = await Promise.all(templateFiles.map(async (f) => {
  //   return Object.assign(new Template(), JSON.parse(await _package.files.readAsText(f)));
  // }));
  // v.append(ui.list(templates));
}


//tags: panel
//name: Hit Triage
//input: dataframe table
//output: widget result
//condition: table.tags.get("HitTriage") != null
export function tableWidget(table: DG.DataFrame): DG.Widget {
  let template = templates.find((t) => t.name == table.getTag('HitTriage'));
  if (template == null)
    template = new Template();
  let root = ui.div([ui.h1(template.name), ui.label(template.query), ui.bigButton('Publish', () => { alert('Gotcha') })]);
  return new DG.Widget(root);
}
