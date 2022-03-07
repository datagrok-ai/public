import { study } from "../clinical-study";
import * as grok from 'datagrok-api/grok';
import * as ui from "datagrok-api/ui";
import * as DG from 'datagrok-api/dg';
import { AE_BROWSER_VIEW_NAME, TIMELINES_VIEW_NAME } from "../constants/view-names-constants";
import * as sdtmCols from "../constants/columns-constants";
import { AE_END_DAY_FIELD, AE_START_DAY_FIELD, AE_TERM_FIELD, VIEWS_CONFIG } from "../views-config";
import { AEBrowserHelper } from "../helpers/ae-browser-helper";
import { ValidationHelper } from "../helpers/validation-helper";
import { c } from "../package";
import { createValidationErrorsDiv } from "./views-validation-utils";
import { updateDivInnerHTML } from "./utils";

export function createAEBrowserHelper(): any {
    const aeBrowserDf = study.domains.ae.clone();
    const aeBrowserHelper = new AEBrowserHelper(aeBrowserDf);
    const timelinesView = grok.shell.view(TIMELINES_VIEW_NAME) as any;
    if (timelinesView) {
      timelinesView.aeBrowserHelper = aeBrowserHelper;
    }
    aeBrowserDf.onCurrentRowChanged.subscribe(() => {
      aeBrowserHelper.currentSubjId = aeBrowserDf.get(sdtmCols.SUBJECT_ID, aeBrowserDf.currentRowIdx);
      aeBrowserHelper.currentAeDay = aeBrowserDf.get(VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_START_DAY_FIELD], aeBrowserDf.currentRowIdx);
      aeBrowserHelper.propertyPanel();
    })
    return { helper: aeBrowserHelper, df: aeBrowserDf };
  }
  
  export function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/' + view.name;
    grok.shell.addView(view);
    return view;
  }
  
  export function createTableView(
    domainsAndColsToCheck: any,
    viewName: string,
    helpUrl: string,
    createViewHelper: (params: any) => any,
    paramsForHelper?: any) {
    let tableView;
    let viewHelper;
    let validator = new ValidationHelper(domainsAndColsToCheck);
    if (validator.validate()) {
      let { helper, df } = createViewHelper(paramsForHelper);
      tableView = DG.TableView.create(df, false);
      viewHelper = helper
    } else {
      tableView = DG.View.create();
      updateDivInnerHTML(tableView.root, createValidationErrorsDiv(validator.missingDomains, validator.missingColumnsInReqDomains, validator.missingColumnsInOptDomains));
    }
    tableView.name = viewName;
    if (helpUrl) {
      tableView.helpUrl = helpUrl;
    }
    return { helper: viewHelper, view: tableView };
  }
  
  
  export function getTableViewsParams() {
    return {
      [AE_BROWSER_VIEW_NAME]: {
        domainsAndColsToCheck: {
          'req_domains': {
            'ae': {
              'req': [
                VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_TERM_FIELD],
                VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_START_DAY_FIELD],
                VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_END_DAY_FIELD]
              ]
            }
          }
        },
        helpUrl: 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/ae_browser.md',
        createViewHelper: createAEBrowserHelper
      }
    }
  }
  