/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { u2 } from '@datagrok-libraries/utils/src/u2';

import { MolTrackDockerService } from './services/moltrack-docker-service';
import { RegistrationView } from './views/registration-tab';
import { MOLTRACK_ENTITY_LEVEL, MOLTRACK_IS_STATIC_FIELD, SAVED_SEARCHES_NODE, Scope, SEARCH_NODE } from './utils/constants';
import { createSearchNode, createSearchView, getSavedSearches, handleSearchURL, loadSearchFields, molTrackSearchFieldsArr } from './views/search';
import { registerAllData, registerAssayData, updateAllMolTrackSchemas } from './utils/registration-utils';
import { batchView, compoundView, createPath, initBulkRegisterView, initRegisterView } from './utils/view-utils';
import { flattened, getCorporateCompoundIdByExactStructure } from './utils/utils';
import { molTrackPropPanel } from './widgets/moltrack-property-panel';
import { PropertySchemaView } from './views/schema-view';

export const _package = new DG.Package();

//tags: init
export async function init(): Promise<void> {
  try {
    await grok.functions.call('Chem:initChemAutostart');
  } catch (e) {}
  // This will be used for the updated docker setup later.
  const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
  const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
  for (const query of queries)
    await (query.prepare()).call();
}

//name: initDB
export async function initDB() {
  await updateAllMolTrackSchemas();
  await registerAssayData();
  await registerAllData();
}

//tags: app
//name: MolTrack
//meta.icon: images/moltrack.png
//input: string path {meta.url: true; optional: true}
//output: view v
//meta.browsePath: Chem
export async function molTrackApp(path: string): Promise<DG.ViewBase> {
  const url = new URL(window.location.href);
  const corporateCompoundId = url.searchParams.get('corporate_compound_id');
  const corporateBatchId = url.searchParams.get('corporate_batch_id');

  const hasPath = !!path;
  const isCompoundPath = hasPath && path.includes('Compound');
  const isRegisterPath = hasPath && path.includes('Register');
  const isBatchPath = hasPath && path.includes('Batch');
  const isSearchPath = hasPath && (path.includes(SEARCH_NODE) || path.includes(SAVED_SEARCHES_NODE));
  const isBulkPath = hasPath && path.includes('Bulk');

  const setPathAndReturn = (view: DG.View) => {
    view.path = path;
    return view;
  };

  try {
    if (isSearchPath)
      return handleSearchURL(path);
  } catch (e: any) {
    grok.shell.error(e);
  }

  if (isBulkPath)
    return setPathAndReturn(initBulkRegisterView(false));

  if (corporateCompoundId)
    return await compoundView(corporateCompoundId);

  if (corporateBatchId)
    return await batchView(corporateBatchId);

  if (isCompoundPath || isRegisterPath)
    return setPathAndReturn(initRegisterView('Compound', false));

  if (isBatchPath)
    return setPathAndReturn(initRegisterView('Batch', false));

  const appHeader = u2.appHeader({
    iconPath: _package.webRoot + '/images/moltrack.png',
    learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/MolTrack/README.md',
    description: '- Chemical compound registration system\n' +
      '- Analyze assay data\n' +
      '- Find contextual information on molecules.\n',
    appTitle: 'MolTrack',
    appSubTitle: 'Track, analyze, and manage chemical data',
    bottomLine: true,
  });

  const getStatisticsWidget = async () => {
    const rows: any[][] = await Promise.all(
      Object.values(Scope).map(async (entity) => {
        try {
          const df = await grok.functions.call('MolTrack:retrieveEntity', { scope: entity, flatten: true });
          const count = df?.rowCount ?? 0;

          return [
            entity,
            ui.link(count.toString(), () => grok.shell.addTableView(df)),
          ];
        } catch (e) {
          grok.shell.error(`Failed to retrieve ${entity}: ${e}`);
          return [entity, 'Error'];
        }
      }),
    );

    return ui.table(rows, (row) => row);
  };

  const viewRoot = ui.divV([appHeader]);
  if (!isSearchPath) {
    viewRoot.append(ui.wait(async () => {
      const statsWidget = await getStatisticsWidget();
      return ui.div(statsWidget);
    }));
  }

  return DG.View.fromRoot(viewRoot);
}

//input: dynamic treeNode
//input: view browseView
export async function molTrackAppTreeBrowser(appNode: DG.TreeViewGroup, browseView: any) {
  function createRegisterNode(label: string, initView: () => void) {
    appNode.getOrCreateGroup('Register').item(label).onSelected.subscribe(initView);
  }

  appNode.getOrCreateGroup('Register').onSelected.subscribe(() => {
    const view = initRegisterView('Compound');
    view.path = createPath('Register');
  });
  createRegisterNode('Compound', () => initRegisterView('Compound'));
  createRegisterNode('Batch', () => initRegisterView('Batch'));
  createRegisterNode('Bulk...', () => initBulkRegisterView());
  createRegisterNode('Schema', async () => {
    const schemaView = new PropertySchemaView();
    await schemaView.init();
    schemaView.show();
  });

  const excludedScopes = [Scope.ASSAY_RUNS, Scope.ASSAY_RESULTS];

  //search section
  Object.values(Scope)
    .filter((scope) => !excludedScopes.includes(scope))
    .forEach((scope) => createSearchNode(appNode, scope));

  //saved searches section
  const savedSearchesNode = appNode.getOrCreateGroup(SAVED_SEARCHES_NODE);
  Object.values(Scope)
    .filter((scope) => !excludedScopes.includes(scope))
    .forEach((scope) => {
      const entityGroup = savedSearchesNode.getOrCreateGroup(`${scope.charAt(0).toUpperCase()}${scope.slice(1)}`);
      const savedSearches = getSavedSearches(scope);
      Object.keys(savedSearches).forEach((savedSearch) => {
        const savedSearchNode = entityGroup.item(savedSearch);
        savedSearchNode.onSelected
          .subscribe(() => createSearchView(savedSearch, scope, JSON.parse(savedSearches[savedSearch]), true));
      });
    });
}

//name: fetchCompoundProperties
//description: Retrieves all properties defined for the 'compound' scope
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//output: string result
export async function fetchCompoundProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchCompoundProperties();
}

//name: fetchBatchProperties
//description: Retrieves all properties defined for the 'batch' scope
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
//output: string result
export async function fetchBatchProperties(): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.fetchBatchProperties();
}

//name: fetchSchema
//description: Retrieves all dynamic fields
//output: string result
export async function fetchSchema(): Promise<string> {
  await MolTrackDockerService.init();
  const res = await MolTrackDockerService.fetchSchema();
  return JSON.stringify(res);
}

//name: fetchDirectSchema
//description: Retrieves all static fields
//output: string result
export async function fetchDirectSchema(): Promise<string> {
  await MolTrackDockerService.init();
  const res = await MolTrackDockerService.fetchDirectSchema();
  return JSON.stringify(res);
}

//name: getCompoundByCorporateId
//input: string corporateValue
//output: object compound
export async function getCompoundByCorporateId(corporateValue: string) {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.getCompoundByCorporateId(corporateValue);
}

//name: getBatchByCorporateId
//input: string corporateValue
//output: object batch
export async function getBatchByCorporateId(corporateValue: string) {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.getBatchByCorporateId(corporateValue);
}

// name: registerMolTrackProperties
// input: string jsonPayload
// description: Registers compound properties in the MolTrack service using the provided JSON payload
// output: string result
export async function registerMolTrackProperties(jsonPayload: string): Promise<string> {
  await MolTrackDockerService.init();
  return MolTrackDockerService.updateSchema(jsonPayload);
}

//name: registerAssays
//input: string assayPayload
//output: string result
export async function registerAssays(assayPayload: string): Promise<string> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.registerAssay(assayPayload);
}

//name: registerBulk
//input: file csvFile
//input: string scope
//input: string mapping
//input: string errorHandling
//output: dataframe result
export async function registerBulk(
  csvFile: DG.FileInfo,
  scope: string,
  mapping: string,
  errorHandling: string,
): Promise<DG.DataFrame> {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.registerBulk(csvFile, scope, mapping, errorHandling);
}

//name: search
//input: string query
//input: string entityEndpoint
//output: dataframe df
export async function search(query: string, entityEndpoint: string) {
  await MolTrackDockerService.init();
  return await MolTrackDockerService.search(JSON.parse(query), entityEndpoint);
}

//name: retrieveEntity
//input: string scope
//output: dataframe result
export async function retrieveEntity(scope: string, flatten: boolean = true): Promise<DG.DataFrame | undefined> {
  await MolTrackDockerService.init();
  const resultJson = await MolTrackDockerService.retrieveEntity(scope);
  await loadSearchFields();
  const dynamicProps = molTrackSearchFieldsArr ?
    molTrackSearchFieldsArr.filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === scope && !it.options[MOLTRACK_IS_STATIC_FIELD]) : [];

  const processedRes = flatten ?
    resultJson.map((item: any) => flattened(item, dynamicProps)) :
    resultJson;
  return DG.DataFrame.fromObjects(processedRes);
}

//name: Databases | MolTrack
//input: string mol {semType: Molecule}
//tags: panel
//output: widget res
export async function getMoltrackPropPanelByStructure(mol: string): Promise<DG.Widget> {
  const corporateCompoundId = await getCorporateCompoundIdByExactStructure(mol) ?? '';
  const retrievedCompound = await getCompoundByCorporateId(corporateCompoundId);
  return molTrackPropPanel(retrievedCompound, mol);
}

//name: Databases | MolTrack
//input: string id {semType: Grok ID}
//tags: panel
//output: widget res
export async function getMoltrackPropPanelById(id: string): Promise<DG.Widget> {
  await MolTrackDockerService.init();
  const compound = await MolTrackDockerService.getCompoundByCorporateId(id);
  return molTrackPropPanel(compound);
}
