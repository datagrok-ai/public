import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from '@datagrok-libraries/compute-utils/model-catalog';
import {getPackage, getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {formatSaved, listMyModelRuns} from './components/RunStorage/run-storage';
import {browseRuns} from './components/RunStorage/run-browser-dialog';

const LAST_RUNS = 5;

/** Model Hub's context panel: the model's details and my last runs of it. */
export class Compute2ModelHandler extends ModelHandler {
  override renderProperties(func: DG.Func) {
    const acc = ui.accordion('ComputeModel');
    acc.context = func;
    const title = ui.div([
      ui.span([this.renderIcon(func), ui.label(func.friendlyName), ui.contextActions(func), ui.star(func.id)])]);
    if (func.description)
      title.append(ui.div([ui.markdown(func.description)], 'model-catalog-description'));
    title.append(ui.tags(func), ui.wait(() => this.renderLaunch(func)));
    acc.addTitle(title);
    acc.addPane('Details', () => ui.wait(() => this.renderDetails(func)), true);
    acc.addPane('My runs', () => ui.wait(() => this.renderRuns(func)), true);
    // named Dev so that DevTools doesn't add its code snippet pane
    acc.addPane('Dev', () => ui.tableFromMap({Name: func.nqName}, true));
    acc.end();
    return acc.root;
  }

  private async renderLaunch(func: DG.Func) {
    const missingGroups = ModelHandler.getMissingGroups(func, await this.awaitUserGroups());
    return missingGroups.length > 0 ? this.makeMandatoryGroupsInfo(missingGroups) :
      ui.div([ui.bigButton('Launch', () => ModelHandler.openModel(func))]);
  }

  private async renderDetails(func: DG.Func) {
    const missingGroups = ModelHandler.getMissingGroups(func, await this.awaitUserGroups());
    const groups = JSON.parse(func.options['mandatoryUserGroups'] ?
      `${func.options['mandatoryUserGroups']}` : '[]') as {name: string}[];
    const versions = JSON.parse(func.options['versions'] ? `${func.options['versions']}` : '[]') as string[];
    const pkg = getPackage(func);
    const details: {[key: string]: string | HTMLElement | undefined} = {
      'Type': func.options['editor'] === 'Compute2:TreeWizardEditor' ? 'Workflow' : 'Model',
      'Package': pkg ? `${pkg.name} ${pkg.version}` : undefined,
      'Location': [func.options['department'], func.options['HL_process'], func.options['process']]
        .filter((part) => !!part).join(' › '),
      'Status': func.options['dev.status'],
      'Required groups': groups.map(({name}) =>
        `${name} ${missingGroups.some((group) => group.name === name) ? '✗' : '✓'}`).join(', '),
      'Versions': versions.length > 1 ?
        ui.link('Choose version...', () => ModelHandler.openVersionSelectDialog(func)) : undefined,
      'Updated': func.updatedOn?.format('YYYY-MM-DD'),
    };
    return ui.tableFromMap(Object.fromEntries(Object.entries(details).filter(([_, value]) => !!value)));
  }

  private async renderRuns(func: DG.Func) {
    const runs = await listMyModelRuns(func);
    const browse = ui.link(runs.length > LAST_RUNS ? `Show all ${runs.length}...` : 'Browse runs...',
      () => browseRuns(func));
    if (runs.length === 0)
      return ui.divV([ui.divText('No saved runs yet'), browse]);
    return ui.divV([
      ui.table(runs.slice(0, LAST_RUNS), (run) => [
        ui.link(getRunName(run), () => grok.functions.call('Compute2:OpenWorkflowRun', {id: run.id})
          .catch((e) => grok.shell.error(`Couldn't open ${getRunName(run)}: ${e?.message ?? e}`))),
        formatSaved(run),
      ]),
      browse,
    ]);
  }
}
