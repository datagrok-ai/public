/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Template} from "./template";
import {TemplateHandler} from "./template_handler";

export let _package = new DG.Package();
export let templates: Array<Template>;

//tags: init
export function init() {
  DG.ObjectHandler.register(new TemplateHandler());
  grok.functions.onAfterRunAction.subscribe((call) => {
    if (call.func !instanceof DG.DataQuery || templates == undefined)
      return;
    let template = templates.find((t) => t.context == call.context);
    console.log(call, call.context.getVariable('template'), template);
    if (template == null)
      return;
    for (let param in call.func.outputs) {
      let output = call.outputs.get(param);
      if (output instanceof DG.DataFrame)
        output.setTag('HitTriage', template.name);
    }
  });
}

//tags: app
//name: Hit Triage
export async function app() {
  let v = grok.shell.newView('Hit Triage Templates');
  let templateFiles = await _package.files.list('templates/', false);
  templates = await Promise.all(templateFiles.map(async (f) => {
    return Object.assign(new Template(), JSON.parse(await _package.files.readAsText(f)));
  }));
  v.append(ui.list(templates));
}


//tags: panel
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
