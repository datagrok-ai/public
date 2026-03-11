/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {getKnimeClient} from './knime-client-factory';
import {buildWorkflowInputForm} from './workflow-input-form';
import {parseKnimeOutputs, knimeTableToDataFrame, knimeSpecDataToDataFrame} from './data-conversion';
import {pollJobUntilComplete} from './utils';
import '../css/knime-link.css';

export * from './package.g';
export const _package = new DG.Package();

let openedView: DG.View | null = null;

export class PackageFunctions {
  @grok.decorators.app({
    browsePath: 'Compute',
    name: 'KNIME',
  })
  static async knimeLinkApp(): Promise<DG.ViewBase> {
    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/images/knime.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/KnimeLink/README.md',
      description: '- Connect to your KNIME Business Hub.\n' +
        '- Browse available deployments.\n' +
        '- Pass tables, parameters, and files as input.\n' +
        '- Execute workflows and view results as DataFrames.\n',
    });

    const view = DG.View.fromRoot(appHeader);
    view.name = 'KNIME';
    return view;
  }

  @grok.decorators.appTreeBrowser({app: 'KNIME'})
  static async knimeLinkAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
    let client;
    try {
      client = getKnimeClient();
    }
    catch (e: any) {
      grok.shell.warning('Set KNIME base URL in package settings (Manage > Plugins > KnimeLink)');
      return;
    }

    const setBreadcrumbs = (path: string[], tree: DG.TreeViewGroup) => {
      const breadcrumbs = ui.breadcrumbs(path);
      breadcrumbs.onPathClick.subscribe(async (value: string[]) => {
        const actualItem = value[value.length - 1];
        if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
          return;
        if (actualItem === 'KNIME')
          tree.currentItem = tree;
      });
      if (grok.shell.v) {
        if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') {
          const homeIcon = ui.iconFA('home', () => {
            grok.shell.v.close();
            grok.shell.v = DG.View.createByType(DG.VIEW_TYPE.HOME);
          });
          breadcrumbs.root.firstElementChild!.replaceWith(homeIcon);
        }
        const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
        if (viewNameRoot) {
          viewNameRoot.textContent = '';
          viewNameRoot.appendChild(breadcrumbs.root);
        }
      }
    };

    const openExecutionView = (name: string, executionId: string) => {
      if (openedView)
        openedView.close();

      openedView = DG.View.create(name);
      openedView.name = name;

      const rootDiv = ui.div([], 'knime-workflow-view');
      const inputDiv = ui.div([], 'knime-input-form');
      const messagesDiv = ui.div([], 'knime-node-messages');
      const resultsDiv = ui.div([], 'knime-results');

      // Try to get workflow input descriptors
      client.getWorkflowInputs(executionId).then((workflowSpec) => {
        const formResult = workflowSpec.inputs.length > 0
          ? buildWorkflowInputForm(workflowSpec.inputs) : null;
        if (formResult)
          inputDiv.appendChild(formResult.root);
        else
          inputDiv.appendChild(ui.info('This workflow does not have any inputs.'));

        const runButton = ui.bigButton('Run Workflow', async () => {
          ui.empty(messagesDiv);
          ui.empty(resultsDiv);
          resultsDiv.appendChild(ui.loader());
          const nodeMessages = new Set<string>();
          try {
            const inputs = formResult ? await formResult.getInputs() : {};
            const isRestService = executionId.startsWith('rest:');
            let result: import('./types').KnimeExecutionResult;
            let jobId = '';

            if (isRestService && !workflowSpec.hasFileOutputs) {
              // Sync execution — fast, but job is discarded (no output resource download)
              result = await client.executeSyncWorkflow(executionId, inputs);
            }
            else {
              // Async execution — required for file outputs (persistent job for resource download)
              jobId = await client.startAsyncJob(executionId, inputs);
              result = await pollJobUntilComplete(client, jobId, {
                onProgress: (state, message) => {
                  ui.empty(resultsDiv);
                  resultsDiv.appendChild(ui.divText(`Status: ${state}`));
                  if (message) {
                    nodeMessages.add(message);
                    ui.empty(messagesDiv);
                    messagesDiv.appendChild(ui.info(
                      [...nodeMessages].map((m) => ui.divText(m)), 'Node Messages'));
                  }
                },
              });
            }
            ui.empty(resultsDiv);
            await displayResults(resultsDiv, client, jobId, result, openedView!);
          }
          catch (e: any) {
            const msg = e?.message ?? String(e);
            ui.empty(resultsDiv);
            resultsDiv.appendChild(ui.div([ui.divText('Error: ' + msg)], 'knime-error'));
            grok.shell.error(msg);
          }
        });

        const buttonDiv = ui.div([runButton], 'knime-run-button');
        rootDiv.appendChild(inputDiv);
        rootDiv.appendChild(buttonDiv);
        rootDiv.appendChild(messagesDiv);
        rootDiv.appendChild(resultsDiv);
      });

      openedView.root.appendChild(rootDiv);
      grok.shell.addPreview(openedView);
      setBreadcrumbs(['Home', 'KNIME', name], treeNode);
    };

    treeNode.items.length = 0;
    try {
      const deployments = await client.listDeployments();
      if (deployments.length === 0) {
        treeNode.item('No deployments found');
        return;
      }
      for (const dep of deployments) {
        const node = treeNode.item(dep.name);
        node.onSelected.subscribe(async () => openExecutionView(dep.name, dep.id));
      }
    }
    catch (e: any) {
      treeNode.item(`Error loading deployments: ${e?.message ?? e}`);
    }
  }

  @grok.decorators.func({name: 'Execute KNIME Workflow'})
  static async executeKnimeWorkflow(
    workflowId: string,
    @grok.decorators.param({options: {nullable: true}}) inputJson?: string,
    @grok.decorators.param({options: {nullable: true}}) inputTable?: DG.DataFrame,
    @grok.decorators.param({options: {nullable: true}}) tableParamName?: string,
  ): Promise<DG.DataFrame> {
    const client = getKnimeClient();
    const inputs: {[key: string]: any} = {};

    if (inputJson) {
      const parsed = JSON.parse(inputJson);
      for (const key of Object.keys(parsed))
        inputs[key] = parsed[key];
    }

    if (inputTable) {
      const {dataFrameToKnimeTable} = await import('./data-conversion');
      inputs[tableParamName || 'table-input'] = dataFrameToKnimeTable(inputTable);
    }

    const isRestService = workflowId.startsWith('rest:');
    const hasFileInputs = Object.values(inputs).some((v) => v instanceof Blob || v instanceof File);
    let result: import('./types').KnimeExecutionResult;
    let jobId = '';

    if (isRestService && !hasFileInputs) {
      result = await client.executeSyncWorkflow(workflowId, inputs);
    }
    else {
      jobId = await client.startAsyncJob(workflowId, inputs);
      result = await pollJobUntilComplete(client, jobId);
    }

    const {tables, variables} = parseKnimeOutputs(result.outputs);
    if (result.outputTables)
      for (const t of result.outputTables) {
        const df = t?.['table-spec'] && t?.['table-data']
        ? knimeSpecDataToDataFrame(t['table-spec'], t['table-data'])
        : Array.isArray(t) ? knimeTableToDataFrame(t) : knimeTableToDataFrame([t]);
        tables.push(df);
      }

    if (result.fetchedResources)
      for (const {name: fileName, resource} of result.fetchedResources) {
        if (resource.json !== undefined && Array.isArray(resource.json)) {
          const df = knimeTableToDataFrame(resource.json, fileName);
          df.name = fileName;
          tables.push(df);
        }
        else if (resource.text !== undefined) {
          const df = tryParseTextAsDataFrame(resource.text, fileName, resource.contentType);
          if (df)
            tables.push(df);
        }
      }

    if (tables.length > 0)
      return tables[0];
    if (Object.keys(variables).length > 0)
      return DG.DataFrame.fromObjects([variables])!;
    return DG.DataFrame.create();
  }
}

async function displayResults(
  container: HTMLElement,
  client: ReturnType<typeof getKnimeClient>,
  jobId: string,
  result: import('./types').KnimeExecutionResult,
  view?: DG.View,
): Promise<void> {
  // Use jobId from sync execution response as fallback
  const effectiveJobId = jobId || result.jobId || '';
  const {tables, variables} = parseKnimeOutputs(result.outputs);

  if (result.outputTables)
    for (const t of result.outputTables) {
      const df = t?.['table-spec'] && t?.['table-data']
        ? knimeSpecDataToDataFrame(t['table-spec'], t['table-data'])
        : Array.isArray(t) ? knimeTableToDataFrame(t) : knimeTableToDataFrame([t]);
      tables.push(df);
    }

  // Process output resources — use pre-fetched resources (from getJobResult) or fetch now
  const fileOutputs: {name: string; blob: Blob; contentType: string}[] = [];
  if (result.fetchedResources)
    for (const {name: fileName, resource} of result.fetchedResources)
      processOutputResource(resource, fileName, tables, fileOutputs, variables);
  else if (result.outputResources && effectiveJobId)
    for (const resourceId of Object.keys(result.outputResources)) {
      try {
        const resource = await client.fetchOutputResource(effectiveJobId, resourceId);
        const fileName = result.outputResources![resourceId] || resourceId;
        processOutputResource(resource, fileName, tables, fileOutputs, variables);
      }
      catch (e: any) {
        console.warn(`Failed to fetch output resource '${resourceId}':`, e);
        variables[resourceId] = result.outputResources[resourceId];
      }
    }

  const tabs: {[name: string]: HTMLElement} = {};

  const uniqueTabName = (base: string): string => {
    if (!(base in tabs)) return base;
    let n = 2;
    while (`${base} (${n})` in tabs) n++;
    return `${base} (${n})`;
  };

  // Each table -> its own tab
  for (let i = 0; i < tables.length; i++) {
    const df = tables[i];
    const tabName = uniqueTabName(df.name || `Table ${i + 1}`);
    const tabContent = ui.div([], {style: {width: '100%', height: '100%'}});
    const grid = df.plot.grid();
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    tabContent.appendChild(grid.root);
    if (view) {
      const addButton = ui.icons.add(() => grok.shell.addTablePreview(df), 'Open as table view');
      tabContent.appendChild(addButton);
    }
    tabs[tabName] = tabContent;
  }

  // Each file -> its own tab
  for (const file of fileOutputs) {
    const url = URL.createObjectURL(file.blob);
    const ext = guessFileExtension(file.contentType, file.name);
    const fileName = file.name.includes('.') ? file.name : `${file.name}${ext}`;
    const link = ui.link(`Download ${fileName}`, () => {
      const a = document.createElement('a');
      a.href = url;
      a.download = fileName;
      a.click();
    }, `Download output file (${file.contentType})`);
    tabs[uniqueTabName(fileName)] = ui.div([ui.iconFA('file-download'), link], 'knime-file-output');
  }

  // Variables -> one tab
  if (Object.keys(variables).length > 0) {
    const varDiv = ui.div([], 'knime-variables');
    const rows = Object.entries(variables).map(([k, v]) =>
      `<tr><td><strong>${escapeHtml(k)}</strong></td><td>${escapeHtml(typeof v === 'object' ? JSON.stringify(v) : String(v))}</td></tr>`
    ).join('');
    varDiv.innerHTML = `<table><tr><th>Variable</th><th>Value</th></tr>${rows}</table>`;
    tabs['Variables'] = varDiv;
  }

  // // Data-app -> one tab
  // if (result.dataAppUrl) {
  //   const linkDiv = ui.div([], 'knime-data-app');
  //   const link = ui.link('Open Data App in KNIME Hub', result.dataAppUrl, 'Opens the interactive data app in a new tab');
  //   (link as HTMLAnchorElement).target = '_blank';
  //   linkDiv.appendChild(link);
  //   tabs['Data App'] = linkDiv;
  // }

  if (Object.keys(tabs).length === 0) {
    container.appendChild(ui.divText('Workflow completed with no output.'));
    return;
  }

  if (Object.keys(tabs).length === 1) {
    container.appendChild(Object.values(tabs)[0]);
    return;
  }

  const tabControl = ui.tabControl(tabs);
  tabControl.root.style.width = '100%';
  tabControl.root.style.height = '100%';
  container.appendChild(tabControl.root);
}

/** Try to parse text content (CSV/TSV) as a DataFrame. */
function tryParseTextAsDataFrame(text: string, name: string, contentType: string): DG.DataFrame | null {
  const isCsv = contentType.includes('csv') || contentType.includes('comma') || name.endsWith('.csv');
  const isTsv = contentType.includes('tab-separated') || name.endsWith('.tsv');
  if (!isCsv && !isTsv)
    return null;
  try {
    const df = DG.DataFrame.fromCsv(text);
    df.name = name;
    return df;
  }
  catch {
    return null;
  }
}

function processOutputResource(
  resource: import('./types').KnimeOutputResource,
  fileName: string,
  tables: DG.DataFrame[],
  fileOutputs: {name: string; blob: Blob; contentType: string}[],
  variables: {[key: string]: any},
): void {
  if (resource.json !== undefined) {
    if (Array.isArray(resource.json)) {
      const df = knimeTableToDataFrame(resource.json, fileName);
      df.name = fileName;
      tables.push(df);
    }
    else
      variables[fileName] = resource.json;
  }
  else if (resource.text !== undefined) {
    const df = tryParseTextAsDataFrame(resource.text, fileName, resource.contentType);
    if (df)
      tables.push(df);
    else
      fileOutputs.push({name: fileName, blob: new Blob([resource.text], {type: resource.contentType}), contentType: resource.contentType});
  }
  else if (resource.blob)
    fileOutputs.push({name: fileName, blob: resource.blob, contentType: resource.contentType});
}

function guessFileExtension(contentType: string, name: string): string {
  if (name.includes('.'))
    return '';
  if (contentType.includes('csv')) return '.csv';
  if (contentType.includes('tab-separated')) return '.tsv';
  if (contentType.includes('json')) return '.json';
  if (contentType.includes('xml')) return '.xml';
  if (contentType.includes('plain')) return '.txt';
  if (contentType.includes('pdf')) return '.pdf';
  if (contentType.includes('png')) return '.png';
  if (contentType.includes('jpeg')) return '.jpg';
  return '';
}

function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
