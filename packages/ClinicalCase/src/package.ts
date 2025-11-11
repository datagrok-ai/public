/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {StudySummaryView} from './views/study-summary-view';
import {TimelinesView} from './views/timelines-view';
import {PatientProfileView} from './views/patient-profile-view';
import {AdverseEventsView} from './views/adverse-events-view';
import {ValidationView} from './views/validation-view';
import {LaboratoryView} from './views/laboratory-view';
import {AERiskAssessmentView} from './views/ae-risk-assessment-view';
import {SurvivalAnalysisView} from './views/survival-analysis-view';
import {BoxPlotsView} from './views/boxplots-view';
import {MatrixesView} from './views/matrixes-view';
import {TimeProfileView} from './views/time-profile-view';
import {STUDY_ID} from './constants/columns-constants';
import {TreeMapView} from './views/tree-map-view';
import {MedicalHistoryView} from './views/medical-history-view';
import {VisitsView} from './views/visits-view';
import {StudyConfigurationView} from './views/study-config-view';
import {ADVERSE_EVENTS_VIEW_NAME, AE_BROWSER_VIEW_NAME, AE_RISK_ASSESSMENT_VIEW_NAME,
  CORRELATIONS_VIEW_NAME, DISTRIBUTIONS_VIEW_NAME, LABORATORY_VIEW_NAME, MEDICAL_HISTORY_VIEW_NAME,
  PATIENT_PROFILE_VIEW_NAME, QUESTIONNAIRES_VIEW_NAME, STUDY_CONFIGURATIN_VIEW_NAME, SUMMARY_VIEW_NAME,
  SURVIVAL_ANALYSIS_VIEW_NAME, TIMELINES_VIEW_NAME, TIME_PROFILE_VIEW_NAME, TREE_MAP_VIEW_NAME,
  VALIDATION_VIEW_NAME, VISITS_VIEW_NAME} from './constants/view-names-constants';
import {createClinCaseTableView, TABLE_VIEWS_META} from './utils/views-creation-utils';
//import {CohortView} from './views/cohort-view';
import {QuestionnaiesView} from './views/questionnaires-view';
import {ClinCaseTableView, ClinStudyConfig} from './utils/types';
import {domainsToValidate, StudyConfigFileName} from './constants/constants';
import {ClinicalStudy, studies} from './clinical-study';
import {ClinicalCaseViewBase} from './model/ClinicalCaseViewBase';
import '../css/clinical-case.css';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {scripts} from './package-api';
import {ClinicalCaseViewsConfig} from './views-config';
import {Subject} from 'rxjs';

export * from './package.g';

export const _package = new DG.Package();

export let validationRulesList = null;

export const VIEWS: {[key: string]: {[key: string]: DG.ViewBase}} = {};
export const studiesViewsConfigs: {[key: string]: ClinicalCaseViewsConfig} = {};

const loadingStudyData: {[key: string]: boolean} = {};

const CLINICAL_CASE_APP_PATH: string = '/apps/ClinicalCase';

let currentOpenedView: DG.ViewBase | null = null;

const studyLoadedSubject = new Subject<{name: string, loaded: boolean}>();

type CurrentStudyAndView = {
  study: string;
  viewName: string;
}

export let existingStudies: {[key: string]: ClinStudyConfig} | null = null;
const validationNodes: {[key: string]: DG.TreeViewNode} = {};

const domains = (studyId: string, exactDomains?: string[]) =>
  (exactDomains ?? Object.keys(studies[studyId].domains)).map((it) => it.toLocaleLowerCase());
export let c: DG.FuncCall;

function getCurrentStudyAndView(path: string): CurrentStudyAndView {
  let currentStudy = '';
  let currentViewName = '';
  if (path && !path.startsWith('/browse')) {
    const pathSegments = path.split('/');
    currentStudy = pathSegments[1];
    currentViewName = pathSegments.length > 2 ? pathSegments[2] : SUMMARY_VIEW_NAME;
  }
  return {study: currentStudy, viewName: currentViewName};
}

async function clinicalCaseAppTB(treeNode: DG.TreeViewGroup,
  currentStudy: string, currentViewName: string) {
  const loaderDiv = ui.div([], {style: {width: '50px', height: '24px', position: 'relative'}});
  loaderDiv.innerHTML = `<div class="grok-loader"><div></div><div></div><div></div><div></div></div>`;
  const loaderItem = treeNode.item(loaderDiv);
  await loadStudiesFromAppData([]);
  openStudy(treeNode, currentStudy, currentViewName);
  loaderItem.remove();
}

export function openStudy(treeNode: DG.TreeViewGroup,
  currentStudy: string, currentViewName: string, files?: DG.FileInfo[]) {
  let initialViewSelected = false;

  const openStudyNode = async (study: ClinStudyConfig, node: DG.TreeViewGroup) => {
    const summaryView = !VIEWS[study.name] ? null : VIEWS[study.name][SUMMARY_VIEW_NAME];
    if (!summaryView) {
      if (loadingStudyData[study.name]) {
        grok.shell.warning(`Loading data for study ${study.name}`);
        return;
      }
    }
    loadView(study, SUMMARY_VIEW_NAME, node);
  };


  const loadView = async (study: ClinStudyConfig, viewName: string, parentNode: DG.TreeViewGroup) => {
    let view = VIEWS[study.name][viewName];
    let helper: any;
    if (!view) {
      if (Object.keys(TABLE_VIEWS_META).includes(viewName)) { //load table view
        const clinCaseTableView = VIEW_CREATE_FUNC[viewName](study.name, viewName) as ClinCaseTableView;
        view = clinCaseTableView.view;
        helper = clinCaseTableView.helper;
      } else { // load view
        if (viewName === SUMMARY_VIEW_NAME) {
          const errorLinkHandler = () => {
            parentNode.expanded = true;
            parentNode.currentItem = validationNodes[study.name];
          };
          view = VIEW_CREATE_FUNC[viewName](study.name, errorLinkHandler) as DG.ViewBase;
        } else
          view = VIEW_CREATE_FUNC[viewName](study.name) as DG.ViewBase;
      }
    }
    currentOpenedView?.close();
    currentOpenedView = grok.shell.addPreview(view);
    if (view.hasOwnProperty('loaded') && !(view as ClinicalCaseViewBase).loaded)
      (view as ClinicalCaseViewBase).load();
    else
      helper?.propertyPanel();
    view.path =
        `browse${CLINICAL_CASE_APP_PATH}/${study.name}/${viewName.replaceAll(' ', '')}`;
  };

  // eslint-disable-next-line no-unused-vars
  for (const [_, study] of Object.entries(existingStudies)) {
    const node = treeNode.getOrCreateGroup(study.friendlyName ?? study.name, null, false);
    if (!studies[study.name]) {
      studies[study.name] = new ClinicalStudy(study.name);
      VIEWS[study.name] = {};
    }

    node.onSelected.subscribe(async (_) => {
      await openStudyNode(study, node);
    });

    node.onNodeExpanding.subscribe(async (_) => {
      if (loadingStudyData[study.name])
        return;
      for (const viewName of Object.keys(VIEW_CREATE_FUNC)) {
        const viewNode = node.item(viewName);
        if (viewName === VALIDATION_VIEW_NAME)
          validationNodes[study.name] = viewNode;
        viewNode.onSelected.subscribe(() => {
          if (loadingStudyData[study.name]) {
            grok.shell.warning(`Loading data for study ${study.name}`);
            treeNode.currentItem = node;
            return;
          }
          loadView(study, viewName, node);
        });
      }
      await initClinicalStudy(study, files);
      const viewItem = node.items
        .find((node) => node.text === (initialViewSelected ? currentViewName : SUMMARY_VIEW_NAME))?.root;
      viewItem?.click();
      initialViewSelected = false;
    });
  }

  if (currentStudy && !Object.keys(existingStudies).includes(currentStudy))
    grok.shell.error(`Study ${currentStudy} doesn't exist`);
  else if (currentStudy) {
    const studyNode = treeNode.getOrCreateGroup(currentStudy);
    if (currentViewName && !Object.keys(VIEW_CREATE_FUNC).includes(currentViewName)) {
      grok.shell.warning(`${currentViewName} view doesn't exist, opening summary view`);
      currentViewName = SUMMARY_VIEW_NAME;
    } else if (!currentViewName)
      currentViewName = SUMMARY_VIEW_NAME;
    initialViewSelected = true;
    if (studyNode.expanded)
      openStudyNode(existingStudies[currentStudy], studyNode);
    else
      studyNode.expanded = true;
  }
}


export async function initClinicalStudy(study: ClinStudyConfig, studyFiles?: DG.FileInfo[]) {
  if (!studies[study.name].initCompleted) {
    try {
      loadingStudyData[study.name] = true;
      const progressBar = DG.TaskBarProgressIndicator.create(`Reading data for study ${study.name}`);
      const dataLoaded = await readClinicalData(study, undefined, studyFiles);
      if (!dataLoaded)
        studyLoadedSubject.next({name: study.name, loaded: false});
      studies[study.name].init();
      progressBar.close();
      grok.shell.info(`Data for study ${study.name} is ready`);
      studyLoadedSubject.next({name: study.name, loaded: true});
    } catch (e: any) {
      studyLoadedSubject.next({name: study.name, loaded: false});
      throw e;
    } finally {
      delete loadingStudyData[study.name];
    }
  }
}


export async function readClinicalData(study: ClinStudyConfig, domainsToDownLoad?: string[],
  studyFilesToRead?: DG.FileInfo[]): Promise<boolean> {
  const pb = DG.TaskBarProgressIndicator.create(`Reading data for ${study.friendlyName ?? study.name}...`);
  try {
    const studyFiles = studyFilesToRead ?? await _package.files.list(`studies/${study.name}`);
    const domainsList = domains(study.name, domainsToDownLoad);
    const removeExtension = (filename: string) => {
      const lastDotIndex = filename.lastIndexOf('.');
      return lastDotIndex === -1 ? filename : filename.substring(0, lastDotIndex);
    };
    for (let i = 0; i < studyFiles.length; i++) {
      const domainNameWithExt = studyFiles[i].fileName.toLowerCase();
      const isXpt = domainNameWithExt.endsWith('.xpt');
      const domainNameWithoutExt = removeExtension(domainNameWithExt);
      pb.update(i/studyFiles.length * 100, `Reading ${domainNameWithExt}...`);
      if (!studies[study.name].domains[domainNameWithoutExt] && domainsList.includes(domainNameWithoutExt)) {
        let df: DG.DataFrame | null = null;
        if (isXpt) {
          //const bytes = await _package.files.(`studies/${study.name}/${domainNameWithExt}`);
          console.log(`*************** read ${domainNameWithExt}`);
          df = await scripts.readSas(studyFiles[i]);
          console.log(`*************** converted ${domainNameWithExt}`);
        } else
          df = await DG.DataFrame.fromCsv(await studyFiles[i].readAsString());
        if (df) {
          df.name = domainNameWithoutExt;
          studies[study.name].domains[domainNameWithoutExt] = df;
        }
      }
      studiesViewsConfigs[study.name] = new ClinicalCaseViewsConfig();
    }
    if (!studies[study.name].domains.dm) {
      grok.shell.error(`No demographic data found for study ${study.name}`);
      return false;
    }
    return true;
  } finally {
    pb.close();
  }
}


export async function loadStudiesFromAppData(deletedCampaigns: string[]): Promise<{[name: string]: ClinStudyConfig}> {
  if (!existingStudies) {
    const studiesFolders = (await _package.files.list(`studies`))
      .filter((f) => deletedCampaigns.indexOf(f.name) === -1);
    const studiesNamesMap: {[name: string]: ClinStudyConfig} = {};
    for (const folder of studiesFolders) {
      try {
        const studyJson: ClinStudyConfig = JSON.parse(await _package.files
          .readAsText(`studies/${folder.name}/${StudyConfigFileName}`));
        studiesNamesMap[studyJson.name] = studyJson;
      } catch (e) {
        continue;
      }
    }
    existingStudies = studiesNamesMap;
  }
  return existingStudies;
}

export class PackageFunctions {
  @grok.decorators.app({
    'browsePath': 'Clinical',
    'name': 'Clinical Case',
  })
  static async clinicalCaseApp(): Promise<DG.ViewBase | void> {
    const appHeader = u2.appHeader({
      iconPath: _package.webRoot + '/img/clin_case_icon.png',
      learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/ClinicalCase/README.md',
      description:
      '-  Visualize and explore your SDTM data\n' +
      '-  Find patterns and trends in your data\n' +
      '-  Explore data on patient specific and trial specific levels\n' +
      '-  Browse through AEs and related data\n' +
      '-  Validate your SDTM data',
    });

    const studiesHeader = ui.h1('Studies');
    const existingStudies = await loadStudiesFromAppData([]);
    const existingStudiesNames = Object.keys(existingStudies);
    const studiesDiv = ui.divV([]);
    const addStudyLink = (studyName: string) => {
      const studyLink = ui.link(studyName, async () => {
        const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
          .getOrCreateGroup('Clinical').getOrCreateGroup('Clinical Case');
        openStudy(clinicalCaseNode, studyName, SUMMARY_VIEW_NAME);
      }, 'Click to run the study');
      studyLink.style.paddingBottom = '10px';
      studiesDiv.append(studyLink);
    };
    for (const studyName of existingStudiesNames)
      addStudyLink(studyName);

    const importStudyDiv = ui.div('', {style: {position: 'relative'}});
    const importStudyButton = ui.bigButton('Import study', () => {
      const dialog = ui.dialog('Import study');
      let studyConfig: ClinStudyConfig | null = null;
      const filesInput = ui.input.file('Files', {
        directoryInput: true,
        onValueChanged: async () => {
          //check if configuration file is present in the folder
          const studyConfigFile = filesInput.directoryFiles.filter((it) => it.name === StudyConfigFileName);
          if (!studyConfigFile.length)
            // eslint-disable-next-line max-len
            grok.shell.error(`Study folder must include study.json configuration file containing study name, e.g. { "name": "study_name" }`);
          try {
            studyConfig = JSON.parse(await studyConfigFile[0].readAsString());
            if (!studyConfig.name)
              grok.shell.error(`study.json must study name, e.g. { "name": "study_name" }`);
          } catch (e) {
            grok.shell.error(`Error reading study.json`);
          }
          dialog.getButton('OK').disabled = studyConfig === null;
        },
      });
      dialog.add(filesInput.root).show()
        .onOK(() => {
          ui.setUpdateIndicator(importStudyDiv, true, `Loading data for study ${studyConfig.name}`);
          existingStudies[studyConfig.name] = studyConfig;
          studies[studyConfig.name] = new ClinicalStudy(studyConfig.name);
          VIEWS[studyConfig.name] = {};
          const sub = studyLoadedSubject.subscribe((data) => {
            if (data.name === studyConfig.name) {
              ui.setUpdateIndicator(importStudyDiv, false);
              sub.unsubscribe();
              if (data.loaded)
                addStudyLink(studyConfig.name);
            }
          });
          //need setTimeout for dialog not to freeze on cliking OK
          setTimeout(() => {
            openStudy(grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps')
              .getOrCreateGroup('Clinical').getOrCreateGroup('Clinical Case'),
            studyConfig.name, SUMMARY_VIEW_NAME, filesInput.directoryFiles);
          }, 10);
        });
    });
    importStudyButton.style.marginLeft = '0px';
    importStudyDiv.append(importStudyButton);
    const view = DG.View.create();
    view.name = 'Clinical Case';
    view.path = CLINICAL_CASE_APP_PATH;
    view.root.append(ui.divV([
      appHeader,
      studiesHeader,
      studiesDiv,
      importStudyDiv,
    ]));
    return view;
  }

  @grok.decorators.func()
  static async clinicalCaseAppTreeBrowser(treeNode: DG.TreeViewGroup) {
    if (!validationRulesList)
      validationRulesList = await grok.data.loadTable(`${_package.webRoot}tables/validation-rules.csv`);
    const url = new URL(window.location.href);
    const currentStudyAndViewPath = url.pathname.includes(`/browse${CLINICAL_CASE_APP_PATH}`) ?
      url.pathname.replace(`/browse${CLINICAL_CASE_APP_PATH}`, ``) : '';
    const studyAndView = getCurrentStudyAndView(currentStudyAndViewPath);
    await clinicalCaseAppTB(treeNode, studyAndView.study, studyAndView.viewName);
  }

  @grok.decorators.folderViewer({
    'name': 'clinicalCaseFolderLauncher',
  })
  static async clinicalCaseFolderLauncher(
    folder: DG.FileInfo,
    files: DG.FileInfo[]):
    Promise<DG.Widget | undefined> {
    if (files.some((f) => f.fileName.toLowerCase() === 'dm.csv')) {
      const res = await grok.dapi.files.readAsText(`${folder.fullPath}/dm.csv`);
      const table = DG.DataFrame.fromCsv(res);
      const studyId = table.columns.names().includes(STUDY_ID) ? table.get(STUDY_ID, 0) : 'undefined';
      return DG.Widget.fromRoot(ui.div([
        ui.panel([
          ui.divText('Folder contains SDTM data'),
          ui.divText(`Study ID: ${studyId}`)]),
        ui.button('Run ClinicalCase', async () => {
          studies[studyId] = new ClinicalStudy(studyId);
          await Promise.all(files.map(async (file) => {
            if (domains(studyId).includes(file.fileName.toLowerCase())) {
              const df = await grok.data.files.openTable(`${folder.fullPath}/${file.fileName.toLowerCase()}`);
              grok.shell.addTableView(df);
            }
          }));
          grok.functions.call('Clinicalcase:clinicalCaseApp');
        }),
      ]));
    }
  }

  @grok.decorators.fileHandler({
    'ext': 'xpt',
  })
  static async xptFileHandler(
    @grok.decorators.param({'type': 'list'}) file: DG.FileInfo): Promise<DG.DataFrame[]> {
    const res: DG.DataFrame = await scripts.readSas(file);
    res.name = file.name;
    return [res];
  }
}

export const VIEW_CREATE_FUNC: {[key: string]: (studyId: string, args?: any) => DG.ViewBase | ClinCaseTableView} = {

  [SUMMARY_VIEW_NAME]:
    (studyId, errorLinkHandler?: () => void) => new StudySummaryView(SUMMARY_VIEW_NAME, studyId, errorLinkHandler),
  [TIMELINES_VIEW_NAME]: (studyId) => new TimelinesView(TIMELINES_VIEW_NAME, studyId),
  [LABORATORY_VIEW_NAME]: (studyId) => new LaboratoryView(LABORATORY_VIEW_NAME, studyId),
  [PATIENT_PROFILE_VIEW_NAME]: (studyId) => new PatientProfileView(PATIENT_PROFILE_VIEW_NAME, studyId),
  [ADVERSE_EVENTS_VIEW_NAME]: (studyId) => new AdverseEventsView(ADVERSE_EVENTS_VIEW_NAME, studyId),
  [AE_RISK_ASSESSMENT_VIEW_NAME]: (studyId) => new AERiskAssessmentView(AE_RISK_ASSESSMENT_VIEW_NAME, studyId),
  [SURVIVAL_ANALYSIS_VIEW_NAME]: (studyId) => new SurvivalAnalysisView(SURVIVAL_ANALYSIS_VIEW_NAME, studyId),
  [DISTRIBUTIONS_VIEW_NAME]: (studyId) => new BoxPlotsView(DISTRIBUTIONS_VIEW_NAME, studyId),
  [CORRELATIONS_VIEW_NAME]: (studyId) => new MatrixesView(CORRELATIONS_VIEW_NAME, studyId),
  [TIME_PROFILE_VIEW_NAME]: (studyId) => new TimeProfileView(TIME_PROFILE_VIEW_NAME, studyId),
  [TREE_MAP_VIEW_NAME]: (studyId) => new TreeMapView(TREE_MAP_VIEW_NAME, studyId),
  [MEDICAL_HISTORY_VIEW_NAME]: (studyId) => new MedicalHistoryView(MEDICAL_HISTORY_VIEW_NAME, studyId),
  [VISITS_VIEW_NAME]: (studyId) => new VisitsView(VISITS_VIEW_NAME, studyId),
  [STUDY_CONFIGURATIN_VIEW_NAME]: (studyId, addView?: boolean) =>
    new StudyConfigurationView(STUDY_CONFIGURATIN_VIEW_NAME, studyId, addView),
  [VALIDATION_VIEW_NAME]: (studyId) => new ValidationView(VALIDATION_VIEW_NAME, studyId),
  [QUESTIONNAIRES_VIEW_NAME]: (studyId) => new QuestionnaiesView(QUESTIONNAIRES_VIEW_NAME, studyId),
  [AE_BROWSER_VIEW_NAME]: (studyId, viewName) => createClinCaseTableView(studyId, viewName),
};


