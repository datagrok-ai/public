/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {MolTrackDockerService} from './services/moltrack-docker-service';
import {excludedScopes, IViewContainer, MOLTRACK_ENTITY_LEVEL, MOLTRACK_IS_STATIC_FIELD, MOLTRACK_NODE, SAVED_SEARCHES_NODE, Scope, SEARCH_NODE} from './utils/constants';
import {createSavedSearchesSatistics, createSearchExpandableNode, createSearchNode, createSearchView, getSavedSearches, handleSearchURL, loadSearchFields, molTrackSearchFieldsArr} from './views/search';
import {registerAllData, registerAssayData, updateAllMolTrackSchemas} from './utils/registration-utils';
import {batchView, compoundView, createPath, getQuickActionsWidget, getStatisticsWidget, initAssayRegisterView, initBulkRegisterView, initRegisterView, initSchemaView} from './utils/view-utils';
import {flattened} from './utils/utils';
import {molTrackPropPanel} from './widgets/moltrack-property-panel';
import {registerSemanticTypes} from './utils/detectors';

export * from './package.g';
export const _package = new DG.Package();

export class PackageFunctions {
  @grok.decorators.init()
  static async init(): Promise<void> {
    try {
      await grok.functions.call('Chem:initChemAutostart');
    } catch (e) {}
    registerSemanticTypes();

    const dbInitialized = await grok.data.query('MolTrack:checkDBInitialized');
    if (dbInitialized.get('db_initialized', 0))
      return;

    const connection = await grok.dapi.connections.filter('name = "moltrack"').first();
    const queries = await grok.dapi.queries.filter(`connection.id = "${connection.id}"`).list();
    for (const query of queries)
      await (query.prepare()).call();
  }

  @grok.decorators.func({name: 'initDB'})
  static async initDB() {
    await updateAllMolTrackSchemas();
    await registerAssayData();
    await registerAllData();
  }

  @grok.decorators.app({
    name: 'MolTrack',
    meta: {browsePath: 'Chem', icon: 'images/moltrack.png'},
  })
  static async molTrackApp(
    @grok.decorators.param({options: {metaUrl: true, optional: true}}) path: string,
  ): Promise<DG.ViewBase> {
    const url = new URL(window.location.href);
    const corporateCompoundId = url.searchParams.get('corporate_compound_id');
    const corporateBatchId = url.searchParams.get('corporate_batch_id');

    const hasPath = !!path;
    const isCompoundPath = hasPath && path.includes('Compound');
    const isRegisterPath = hasPath && path.includes('Register');
    const isBatchPath = hasPath && path.includes('Batch');
    const isSchemaPath = hasPath && path.includes('Schema');
    const isSearchPath = hasPath && (path.includes(SEARCH_NODE) || path.includes(SAVED_SEARCHES_NODE));
    const isBulkPath = hasPath && path.includes('Bulk');
    const isAssayPath = hasPath && path.includes('Assay');

    const setPathAndReturn = (registerView: IViewContainer): DG.ViewBase => {
      const {view} = registerView;
      view.path = path;
      return view;
    };

    try {
      if (isSearchPath)
        return await handleSearchURL(path);
    } catch (e: any) {
      grok.shell.error(e);
    }

    if (isBulkPath)
      return setPathAndReturn(initBulkRegisterView());

    if (isAssayPath)
      return setPathAndReturn(initAssayRegisterView());

    if (corporateCompoundId)
      return await compoundView(corporateCompoundId);

    if (corporateBatchId)
      return await batchView(corporateBatchId);

    if (isCompoundPath || isRegisterPath)
      return setPathAndReturn(initRegisterView('Compound', false));

    if (isBatchPath)
      return setPathAndReturn(initRegisterView('Batch', false));

    if (isSchemaPath)
      return setPathAndReturn(await initSchemaView());

    const statisticsWidget = await getStatisticsWidget(createSearchView);
    const quickActionsWidget = getQuickActionsWidget();
    const viewRoot = ui.divH([statisticsWidget, quickActionsWidget], {style: {gap: '50px'}});
    const view = await createSearchExpandableNode([], async () => viewRoot, false);
    view.name = MOLTRACK_NODE;
    return view;
  }

  @grok.decorators.appTreeBrowser({app: 'MolTrack'})
  static async molTrackAppTreeBrowser(
    @grok.decorators.param({type: 'dynamic'}) appNode: DG.TreeViewGroup,
    @grok.decorators.param({type: 'view'}) browseView: any,
  ) {
    function createRegisterNode(label: string, initView: () => IViewContainer | Promise<IViewContainer>) {
      appNode.getOrCreateGroup('Register').item(label).onSelected.subscribe(async () => {
        const registerView = await initView();
        registerView.show?.();
        return registerView.view;
      });
    }

    appNode.getOrCreateGroup('Register').onSelected.subscribe(() => {
      const view = initRegisterView('Compound');
      view.path = createPath('Register');
    });
    createRegisterNode('Compound', () => initRegisterView('Compound'));
    createRegisterNode('Batch', () => initRegisterView('Batch'));
    createRegisterNode('Assay', () => initAssayRegisterView());
    createRegisterNode('Bulk', () => initBulkRegisterView());
    createRegisterNode('Schema', async () => await initSchemaView());

    const searchNode = appNode.getOrCreateGroup(SEARCH_NODE);
    const searchableScopes = Object.values(Scope)
      .filter((scope) => !excludedScopes.includes(scope));
    searchNode.onSelected.subscribe(() =>
      createSearchExpandableNode([SEARCH_NODE], () => getStatisticsWidget(createSearchView)));

    //search section
    searchableScopes.forEach((scope) => createSearchNode(appNode, scope));

    //saved searches section
    const savedSearchesNode = appNode.getOrCreateGroup(SAVED_SEARCHES_NODE);
    savedSearchesNode.onSelected.subscribe(() => createSearchExpandableNode([SAVED_SEARCHES_NODE], () => createSavedSearchesSatistics(undefined)));
    Object.values(Scope)
      .filter((scope) => !excludedScopes.includes(scope))
      .forEach((scope) => {
        const entityGroup = savedSearchesNode.getOrCreateGroup(`${scope.charAt(0).toUpperCase()}${scope.slice(1)}`);
        entityGroup.onSelected.subscribe(() => createSearchExpandableNode([SAVED_SEARCHES_NODE, scope], () => createSavedSearchesSatistics(scope)));
        const savedSearches = getSavedSearches(scope);
        Object.keys(savedSearches).forEach((savedSearch) => {
          const savedSearchNode = entityGroup.item(savedSearch);
          savedSearchNode.onSelected
            .subscribe(() => createSearchView(savedSearch, scope, JSON.parse(savedSearches[savedSearch]), true));
        });
      });
  }

  @grok.decorators.func({
    name: 'fetchCompoundProperties',
    description: 'Retrieves all properties defined for the \'compound\' scope',
    cache: 'all',
    cacheInvalidateOn: '0 0 1 * *',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async fetchCompoundProperties(): Promise<string> {
    return await MolTrackDockerService.fetchCompoundProperties();
  }

  @grok.decorators.func({
    name: 'fetchBatchProperties',
    description: 'Retrieves all properties defined for the \'batch\' scope',
    cache: 'all',
    cacheInvalidateOn: '0 0 1 * *',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async fetchBatchProperties(): Promise<string> {
    return await MolTrackDockerService.fetchBatchProperties();
  }

  @grok.decorators.func({
    name: 'fetchSchema',
    description: 'Retrieves all dynamic fields',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async fetchSchema(): Promise<string> {
    const res = await MolTrackDockerService.fetchSchema();
    return JSON.stringify(res);
  }

  @grok.decorators.func({
    name: 'fetchDirectSchema',
    description: 'Retrieves all static fields',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async fetchDirectSchema(): Promise<string> {
    const res = await MolTrackDockerService.fetchDirectSchema();
    return JSON.stringify(res);
  }

  @grok.decorators.func({
    name: 'getCompoundByCorporateId',
    outputs: [{type: 'object', name: 'compound'}],
  })
  static async getCompoundByCorporateId(
    @grok.decorators.param({}) corporateValue: string,
  ) {
    return await MolTrackDockerService.getCompoundByCorporateId(corporateValue);
  }

  @grok.decorators.func({
    name: 'getBatchByCorporateId',
    outputs: [{type: 'object', name: 'batch'}],
  })
  static async getBatchByCorporateId(
    @grok.decorators.param({}) corporateValue: string,
  ) {
    return await MolTrackDockerService.getBatchByCorporateId(corporateValue);
  }

  @grok.decorators.func({
    name: 'registerMolTrackProperties',
    description: 'Registers compound properties in the MolTrack service using the provided JSON payload',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async registerMolTrackProperties(
    @grok.decorators.param({}) jsonPayload: string,
  ): Promise<string> {
    return MolTrackDockerService.updateSchema(jsonPayload);
  }

  @grok.decorators.func({
    name: 'registerAssays',
    outputs: [{type: 'string', name: 'result'}],
  })
  static async registerAssays(
    @grok.decorators.param({}) assayPayload: string,
  ): Promise<string> {
    return await MolTrackDockerService.registerAssay(assayPayload);
  }

  @grok.decorators.func({
    name: 'registerBulk',
    outputs: [{type: 'dataframe', name: 'result'}],
  })
  static async registerBulk(
    @grok.decorators.param({type: 'file'}) csvFile: DG.FileInfo,
    @grok.decorators.param({}) scope: string,
    @grok.decorators.param({}) mapping: string,
    @grok.decorators.param({}) errorHandling: string,
  ): Promise<DG.DataFrame> {
    return await MolTrackDockerService.registerBulk(csvFile, scope, mapping, errorHandling);
  }

  @grok.decorators.func({
    name: 'search',
    outputs: [{type: 'dataframe', name: 'df'}],
  })
  static async search(
    @grok.decorators.param({}) query: string,
    @grok.decorators.param({}) entityEndpoint: string,
  ) {
    return await MolTrackDockerService.search(JSON.parse(query), entityEndpoint);
  }

  @grok.decorators.func({
    name: 'advancedSearch',
    description: 'Performs a structured MolTrack compound search. The caller must provide the "output": a list of fully-qualified field names to return; "filter": a structured filter tree describing search conditions.',
    outputs: [{type: 'dataframe', name: 'searchResult'}],
  })
  static async advancedSearch(
    @grok.decorators.param({type: 'list', options: {description: 'List of fields to return. All fields must use MolTrack notation. Valid formats: "<table>.<field>" (direct database column) or "<table>.details.<property>" (dynamic detail property). Valid direct compound fields include (non-exhaustive): "compounds.created_at","compounds.updated_at","compounds.created_by","compounds.updated_by","compounds.canonical_smiles","compounds.original_molfile","compounds.molregno","compounds.inchi","compounds.inchikey","compounds.formula","compounds.hash_mol","compounds.hash_tautomer","compounds.hash_canonical_smiles","compounds.hash_no_stereo_smiles","compounds.hash_no_stereo_tautomer","compounds.is_archived".'}}) outputFields: string[],
    @grok.decorators.param({type: 'object', options: {description: 'A structured boolean filter supporting simple and nested conditions. Simple format: {"field":"<field_name>","operator":"<operator>","value":<value>,"threshold":<number|null>} where threshold is required only for "IS SIMILAR". String operators: "=","!=","IN","STARTS WITH","ENDS WITH","LIKE","CONTAINS". Numeric operators: "<",">","<=",">=","RANGE" (expects {value:[min,max]}). Datetime operators: "BEFORE","AFTER","ON" (ISO 8601). Molecular operators (only for compounds.structure): "IS SIMILAR" (requires SMILES + numeric similarity threshold), "IS SUBSTRUCTURE OF", "HAS SUBSTRUCTURE". Complex nested conditions: {"operator":"AND"|"OR","conditions":[...]}. Notes: nested conditions may be arbitrarily deep, each branch must be simple or complex. Examples: simple {"field":"compounds.formula","operator":"=","value":"C6H6"}, complex {"operator":"AND","conditions":[{"field":"compounds.created_at","operator":"AFTER","value":"2025-01-01"},{"operator":"OR","conditions":[{"field":"compounds.formula","operator":"=","value":"C6H6"},{"field":"compounds.structure","operator":"IS SIMILAR","value":"c1ccccc1","threshold":0.8}]}]}.'}}) filter: any,
  ): Promise<DG.DataFrame> {
    const level = 'compounds';
    const molTrackSearchQuery = {
      level: level,
      output: outputFields,
      filter: filter,
      output_format: 'json',
    };
    const searchResultJson = await MolTrackDockerService.search(molTrackSearchQuery, level);
    return DG.DataFrame.fromObjects(searchResultJson.data)!;
  }

  @grok.decorators.func({
    name: 'retrieveEntity',
    outputs: [{type: 'dataframe', name: 'result'}],
  })
  static async retrieveEntity(
    @grok.decorators.param({}) scope: string,
  ): Promise<DG.DataFrame | undefined> {
    const flatten: boolean = true;
    const resultJson = await MolTrackDockerService.retrieveEntity(scope);
    await loadSearchFields();
    const dynamicProps = molTrackSearchFieldsArr ?
      molTrackSearchFieldsArr.filter((it) => it.options[MOLTRACK_ENTITY_LEVEL] === scope && !it.options[MOLTRACK_IS_STATIC_FIELD]) : [];

    const processedRes = flatten ?
      resultJson.map((item: any) => flattened(item, dynamicProps)) :
      resultJson;
    return DG.DataFrame.fromObjects(processedRes);
  }

  @grok.decorators.panel({
    name: 'Databases | MolTrack',
    outputs: [{type: 'widget', name: 'res'}],
  })
  static async getMoltrackPropPanelById(
    @grok.decorators.param({type: 'semantic_value', options: {semType: 'Grok ID'}}) id: DG.SemanticValue,
  ): Promise<DG.Widget> {
    const {value: idValue} = id;
    const compound = await MolTrackDockerService.getCompoundByCorporateId(idValue);
    return molTrackPropPanel(compound, id.cell.column);
  }
}
