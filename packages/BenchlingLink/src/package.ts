/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { queryAASequences, postAASequence } from './aaSequencesApi';
import { queryDNASequences, postDNASequence } from './dnaSequencesApi';
import { queryAssayResults, queryAssayRuns, postAssayResult, postAssayRun } from './assayApi';
import { queryMolecules, postMolecule } from './moleculesApi';
import { dataFrameFromObjects } from './utils';
import {u2} from "@datagrok-libraries/utils/src/u2";
import { queryPlates, postPlate } from './platesApi';
import { queryMixtures, postMixture } from './mixturesApi';
import { queryDnaOligos, postDnaOligo } from './dnaOligosApi';

export * from './package.g';
export const _package = new DG.Package();

const STORAGE_NAME = 'BenchlingLinkFuncEditor';
let openedView: DG.View | null = null;

export class PackageFunctions {
  @grok.decorators.app({
    browsePath: 'Chem',
    name: 'Benchling'
  })
  static async benchlingLinkApp(): Promise<DG.ViewBase> {

    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/images/benchling.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/BenchlingLink/README.md',
      description: '- Integrate with your Benchling account.\n' +
        '- Analyze assay data.\n' +
        '- Find contextual information on molecules and sequences.\n' +
        '- Create entities and post them to the registry.\n' +
        '- Browse the tenant content.\n'
    });

    const view = DG.View.fromRoot(appHeader);
    view.name = 'Benchling';
    return view;
  }

  @grok.decorators.appTreeBrowser({app: 'Benchling'})
  static async benchlingLinkAppTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
    function createFuncEditorView(funcName: string, v: DG.View) {
      let func = DG.Func.byName(funcName);
      let editorDiv = ui.div();
      let gridDiv = ui.div();
      let root = ui.splitV([editorDiv, gridDiv], {style: {height: '100%', width: '100%'}}, true);
      gridDiv.setAttribute('style', 'overflow:hidden !important');
      let df: DG.DataFrame | null;

      const addToWorkspaceButton = ui.icons.add(() => {
        if (df) {
          grok.shell.addTablePreview(df);
        }
      }, 'Add results to workspace');
      v.setRibbonPanels([[addToWorkspaceButton]]);

      // Restore state if exists
      const saved = grok.userSettings.getValue(STORAGE_NAME, funcName);
      const funcCall = func.prepare(saved ? JSON.parse(saved) : null);
      renderEditorAndRun();

      async function renderEditorAndRun() {
        ui.empty(editorDiv);
        let form = await DG.InputForm.forFuncCall(funcCall);
        form.root.style.flexWrap = 'wrap';
        form.root.style.height = '100%';
        form.root.style.maxWidth = 'unset';
        DG.debounce(form.onInputChanged, 1000).subscribe(() => {
          runFunc();
        })
        editorDiv.appendChild(form.root);
        runFunc();
      }

      async function runFunc() {
        //save request settings
        const params: {[key: string]: any} = {};
        Object.keys(funcCall.inputParams).forEach((paramName) => {
          if (funcCall.inputParams[paramName].value)
            params[paramName] = funcCall.inputParams[paramName].value;
        });
        grok.userSettings.add(STORAGE_NAME, funcName, JSON.stringify(params));

        //run function
        ui.empty(gridDiv);
        try {
          df = (await funcCall.call()).getOutputParamValue();
          let grid = df!.plot.grid();
          grid.root.style.width = '100%';
          grid.root.style.height = '95%';
          gridDiv.appendChild(grid.root);
        } catch (e) {
          gridDiv.appendChild(ui.divText('Error: ' + e));
        }
      }

      return root;
    }

    const addBenchlingView = (viewName: string, funcName: string) => {
      if (openedView)
        openedView.close();
      openedView = DG.View.create(viewName);
      openedView.name = viewName;
      openedView.root.appendChild(createFuncEditorView(funcName, openedView));
      const path = ['Home', 'Benchling', viewName];
      grok.shell.addPreview(openedView);
      setBreadcrumbs(path, treeNode);
    }

    const setBreadcrumbs = (path: string[], tree: DG.TreeViewGroup) => {
      const breadcrumbs = ui.breadcrumbs(path);
      breadcrumbs.onPathClick.subscribe(async (value) => {
        const actualItem = value[value.length - 1];
        if (actualItem === breadcrumbs.path[breadcrumbs.path.length - 1])
          return;
        if (actualItem === 'Benchling')
          tree.currentItem = tree;
      });
      if (grok.shell.v) {
        if (breadcrumbs.path.length !== 0 && breadcrumbs.path[0] === 'Home') { // integrate it to the actual breadcrumbs element
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
    }

    try {
      treeNode.items.length = 0;
      const aaSequenceNode = treeNode.item('AA Sequences');
      aaSequenceNode.onSelected.subscribe(async () => {
        addBenchlingView('AA Sequences', 'BenchlingLink:getAASequences');
      });

      const dnaSequencesNode = treeNode.item('DNA Sequences');
      dnaSequencesNode.onSelected.subscribe(async () => {
        addBenchlingView('DNA Sequences', 'BenchlingLink:getDNASequences');
      });

      const assayResultsNode = treeNode.item('Assay Results');
      assayResultsNode.onSelected.subscribe(async () => {
        addBenchlingView('Assay Results', 'BenchlingLink:getAssayResults');
      });

      const assayRunsNode = treeNode.item('Assay Runs');
      assayRunsNode.onSelected.subscribe(async () => {
        addBenchlingView('Assay Runs', 'BenchlingLink:getAssayRuns');
      });

      const moleculesNode = treeNode.item('Molecules');
      moleculesNode.onSelected.subscribe(async () => {
        addBenchlingView('Molecules', 'BenchlingLink:getMolecules');
      });

      const platesNode = treeNode.item('Plates');
      platesNode.onSelected.subscribe(async () => {
        addBenchlingView('Plates', 'BenchlingLink:getPlates');
      });

      const mixturesNode = treeNode.item('Mixtures');
      mixturesNode.onSelected.subscribe(async () => {
        addBenchlingView('Mixtures', 'BenchlingLink:getMixtures');
      });

      const dnaOligosNode = treeNode.item('DNA Oligos');
      dnaOligosNode.onSelected.subscribe(async () => {
        addBenchlingView('DNA Oligos', 'BenchlingLink:getDnaOligos');
      });

      const projectsNode = treeNode.item('Projects');
      projectsNode.onSelected.subscribe(async () => {
        addBenchlingView('Projects', 'BenchlingLink:getProjects');
      });
    } catch (e: any) {
      grok.shell.error(e?.message ?? e);
    }
  }

  @grok.decorators.func({name: 'Get AA Sequences'})
  static async getAASequences(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) aminoAcids?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      aminoAcids,
      folderId,
      mentionedIn,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions,
      ids,
      entityRegistryIds_anyOf,
      names_anyOf,
      names_anyOf_caseSensitive,
      creatorIds,
      authorIds_anyOf,
      returning,
    };
    const aaSequences = await queryAASequences(params);
    return aaSequences;
  }

  @grok.decorators.func({name: 'Get DNA Sequences'})
  static async getDNASequences(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) bases?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      bases,
      folderId,
      mentionedIn,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions,
      ids,
      entityRegistryIds_anyOf,
      names_anyOf,
      names_anyOf_caseSensitive,
      creatorIds,
      authorIds_anyOf,
      returning,
    };
    const dnaSequences = await queryDNASequences(params);
    return dnaSequences;
  }

  @grok.decorators.func({name: 'Get Assay Results'})
  static async getAssayResults(
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_lt?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_gt?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_lte?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt_gte?: string,
    @grok.decorators.param({options: {nullable: true}}) minCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) maxCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) entityIds?: string,
    @grok.decorators.param({options: {nullable: true}}) storageIds?: string,
    @grok.decorators.param({options: {nullable: true}}) assayRunIds?: string,
    @grok.decorators.param({options: {nullable: true}}) automationOutputProcessorId?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_lt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_gt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_lte?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt_gte?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string
  ): Promise<DG.DataFrame> {
    const params = {
      schemaId,
      createdAt_lt,
      createdAt_gt,
      createdAt_lte,
      createdAt_gte,
      minCreatedTime,
      maxCreatedTime,
      sort,
      entityIds,
      storageIds,
      assayRunIds,
      automationOutputProcessorId,
      ids,
      modifiedAt_lt,
      modifiedAt_gt,
      modifiedAt_lte,
      modifiedAt_gte,
      archiveReason,
    };
    const assayResults = await queryAssayResults(params);
    return assayResults;
  }

  @grok.decorators.func({name: 'Get Assay Runs'})
  static async getAssayRuns(
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) minCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) maxCreatedTime?: number,
    @grok.decorators.param({options: {nullable: true}}) ids?: string
  ): Promise<DG.DataFrame> {
    const params = {
      schemaId,
      minCreatedTime,
      maxCreatedTime,
      ids,
    };
    const assayRuns = await queryAssayRuns(params);
    return assayRuns;
  }

  @grok.decorators.func({name: 'Create AA Sequence'})
  static async createAASequence(
    name: string,
    aminoAcids: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
  ): Promise<DG.DataFrame> {
    // Parse JSON string inputs for array/object fields if provided
    const body: any = { name, aminoAcids };
    if (aliases) body.aliases = JSON.parse(aliases);
    if (annotations) body.annotations = JSON.parse(annotations);
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    if (fields) body.fields = JSON.parse(fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    const result = await postAASequence(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Create DNA Sequence'})
  static async createDNASequence(
    name: string,
    bases: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, bases };
    if (aliases) body.aliases = JSON.parse(aliases);
    if (annotations) body.annotations = JSON.parse(annotations);
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    if (fields) body.fields = JSON.parse(fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    const result = await postDNASequence(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Create Assay Result'})
  static async createAssayResult(
    schemaId: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) entityIds?: string,
    @grok.decorators.param({options: {nullable: true}}) storageIds?: string,
    @grok.decorators.param({options: {nullable: true}}) assayRunId?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { schemaId };
    if (fields) body.fields = JSON.parse(fields);
    if (entityIds) body.entityIds = JSON.parse(entityIds);
    if (storageIds) body.storageIds = JSON.parse(storageIds);
    if (assayRunId) body.assayRunId = assayRunId;
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    const result = await postAssayResult(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Create Assay Run'})
  static async createAssayRun(
    schemaId: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { schemaId };
    if (fields) body.fields = JSON.parse(fields);
    if (name) body.name = name;
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    const result = await postAssayRun(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Get Molecules'})
  static async getMolecules(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) chemicalSubstructure_mol?: string,
    @grok.decorators.param({options: {nullable: true}}) chemicalSubstructure_smiles?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      folderId,
      mentionedIn,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions,
      ids,
      entityRegistryIds_anyOf,
      names_anyOf,
      authorIds_anyOf,
      chemicalSubstructure_mol,
      chemicalSubstructure_smiles,
    };
    return await queryMolecules(params);
  }

  @grok.decorators.func({name: 'Create Molecule'})
  static async createMolecule(
    name: string,
    smiles: string,
    @grok.decorators.param({options: {nullable: true}}) formula?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, smiles };
    if (formula) body.formula = formula;
    const result = await postMolecule(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Get Projects'})
  static async getProjects(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      archiveReason,
      ids,
      name,
    };
    const projects = await (await import('./projectsApi')).queryProjects(params);
    return projects;
  }

  @grok.decorators.func({name: 'Get Plates'})
  static async getPlates(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_gte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_gt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_lte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyPositions_lt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_gte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_gt?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_lte?: number,
    @grok.decorators.param({options: {nullable: true}}) emptyContainers_lt?: number,
    @grok.decorators.param({options: {nullable: true}}) ancestorStorageId?: string,
    @grok.decorators.param({options: {nullable: true}}) storageContentsId?: string,
    @grok.decorators.param({options: {nullable: true}}) storageContentsIds?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) barcodes?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      schemaId,
      schemaFields,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      emptyPositions,
      emptyPositions_gte,
      emptyPositions_gt,
      emptyPositions_lte,
      emptyPositions_lt,
      emptyContainers,
      emptyContainers_gte,
      emptyContainers_gt,
      emptyContainers_lte,
      emptyContainers_lt,
      ancestorStorageId,
      storageContentsId,
      storageContentsIds,
      archiveReason,
      ids,
      barcodes,
      names_anyOf,
      names_anyOf_caseSensitive,
      returning,
      creatorIds,
    };
    return await queryPlates(params);
  }

  @grok.decorators.func({name: 'Create Plate'})
  static async createPlate(
    name: string,
    schemaId: string,
    @grok.decorators.param({options: {nullable: true}}) barcode?: string,
    @grok.decorators.param({options: {nullable: true}}) containerSchemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) parentStorageId?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) wells?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, schemaId };
    if (barcode) body.barcode = barcode;
    if (containerSchemaId) body.containerSchemaId = containerSchemaId;
    if (fields) body.fields = JSON.parse(fields);
    if (parentStorageId) body.parentStorageId = parentStorageId;
    if (projectId) body.projectId = projectId;
    if (wells) body.wells = JSON.parse(wells);
    const result = await postPlate(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Get Mixtures'})
  static async getMixtures(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) ingredientComponentEntityIds?: string,
    @grok.decorators.param({options: {nullable: true}}) ingredientComponentEntityIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      folderId,
      mentionedIn,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions,
      ids,
      names_anyOf,
      names_anyOf_caseSensitive,
      entityRegistryIds_anyOf,
      ingredientComponentEntityIds,
      ingredientComponentEntityIds_anyOf,
      authorIds_anyOf,
    };
    return await queryMixtures(params);
  }

  @grok.decorators.func({name: 'Create Mixture'})
  static async createMixture(
    name: string,
    ingredients: string,
    schemaId: string,
    units: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) amount?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryId?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, schemaId, units };
    body.ingredients = JSON.parse(ingredients);
    if (aliases) body.aliases = JSON.parse(aliases);
    if (amount) body.amount = amount;
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    if (entityRegistryId) body.entityRegistryId = entityRegistryId;
    if (fields) body.fields = JSON.parse(fields);
    if (folderId) body.folderId = folderId;
    const result = await postMixture(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }

  @grok.decorators.func({name: 'Get DNA Oligos'})
  static async getDnaOligos(
    @grok.decorators.param({options: {nullable: true}}) sort?: string,
    @grok.decorators.param({options: {nullable: true}}) createdAt?: string,
    @grok.decorators.param({options: {nullable: true}}) modifiedAt?: string,
    @grok.decorators.param({options: {nullable: true}}) name?: string,
    @grok.decorators.param({options: {nullable: true}}) nameIncludes?: string,
    @grok.decorators.param({options: {nullable: true}}) bases?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) mentionedIn?: string,
    @grok.decorators.param({options: {nullable: true}}) projectId?: string,
    @grok.decorators.param({options: {nullable: true}}) registryId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaFields?: string,
    @grok.decorators.param({options: {nullable: true}}) archiveReason?: string,
    @grok.decorators.param({options: {nullable: true}}) mentions?: string,
    @grok.decorators.param({options: {nullable: true}}) ids?: string,
    @grok.decorators.param({options: {nullable: true}}) entityRegistryIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) names_anyOf_caseSensitive?: string,
    @grok.decorators.param({options: {nullable: true}}) creatorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds_anyOf?: string,
    @grok.decorators.param({options: {nullable: true}}) returning?: string,
    @grok.decorators.param({options: {nullable: true}}) customNotationId?: string,
  ): Promise<DG.DataFrame> {
    const params = {
      sort,
      createdAt,
      modifiedAt,
      name,
      nameIncludes,
      bases,
      folderId,
      mentionedIn,
      projectId,
      registryId,
      schemaId,
      schemaFields,
      archiveReason,
      mentions,
      ids,
      entityRegistryIds_anyOf,
      names_anyOf,
      names_anyOf_caseSensitive,
      creatorIds,
      authorIds_anyOf,
      returning,
      customNotationId,
    };
    return await queryDnaOligos(params);
  }

  @grok.decorators.func({name: 'Create DNA Oligo'})
  static async createDnaOligo(
    name: string,
    bases: string,
    @grok.decorators.param({options: {nullable: true}}) aliases?: string,
    @grok.decorators.param({options: {nullable: true}}) annotations?: string,
    @grok.decorators.param({options: {nullable: true}}) authorIds?: string,
    @grok.decorators.param({options: {nullable: true}}) customFields?: string,
    @grok.decorators.param({options: {nullable: true}}) fields?: string,
    @grok.decorators.param({options: {nullable: true}}) folderId?: string,
    @grok.decorators.param({options: {nullable: true}}) schemaId?: string,
    @grok.decorators.param({options: {nullable: true}}) helm?: string,
  ): Promise<DG.DataFrame> {
    const body: any = { name, bases };
    if (aliases) body.aliases = JSON.parse(aliases);
    if (annotations) body.annotations = JSON.parse(annotations);
    if (authorIds) body.authorIds = JSON.parse(authorIds);
    if (customFields) body.customFields = JSON.parse(customFields);
    if (fields) body.fields = JSON.parse(fields);
    if (folderId) body.folderId = folderId;
    if (schemaId) body.schemaId = schemaId;
    if (helm) body.helm = helm;
    const result = await postDnaOligo(body);
    return dataFrameFromObjects([result]) ?? DG.DataFrame.create();
  }
}
