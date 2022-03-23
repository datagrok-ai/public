import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import { ClinicalCaseViewBase } from "../model/ClinicalCaseViewBase";
import { AE_END_DAY_FIELD, AE_START_DAY_FIELD, AE_TERM_FIELD, CON_MED_END_DAY_FIELD, CON_MED_NAME_FIELD, CON_MED_START_DAY_FIELD, INV_DRUG_END_DAY_FIELD, INV_DRUG_NAME_FIELD, INV_DRUG_START_DAY_FIELD, MH_TERM_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from "../views-config";
import { ACT_TRT_ARM, AE_DECOD_TERM, AE_END_DAY, AE_END_DAY_CALCULATED, AE_START_DAY, AE_START_DAY_CALCULATED, AE_TERM, CON_MED_DECOD, CON_MED_END_DAY, CON_MED_END_DAY_CALCULATED, CON_MED_START_DAY, CON_MED_START_DAY_CALCULATED, CON_MED_TRT, INV_DRUG_END_DAY, INV_DRUG_END_DAY_CALCULATED, INV_DRUG_ID, INV_DRUG_NAME, INV_DRUG_START_DAY, INV_DRUG_START_DAY_CALCULATED, MH_DECOD_TERM, MH_TERM, PLANNED_TRT_ARM } from "../constants/columns-constants";
import { updateDivInnerHTML } from "../utils/utils";
import { VIEWS } from "../constants/constants";
import { addView, createTableView, getTableViewsParams } from "../utils/views-creation-utils";

export class StudyConfigurationView extends ClinicalCaseViewBase {

    choices = {
        [TRT_ARM_FIELD]: [ACT_TRT_ARM, PLANNED_TRT_ARM],
        [AE_TERM_FIELD]: [AE_TERM, AE_DECOD_TERM],
        [AE_START_DAY_FIELD]: [AE_START_DAY, AE_START_DAY_CALCULATED],
        [AE_END_DAY_FIELD]: [AE_END_DAY, AE_END_DAY_CALCULATED],
        [INV_DRUG_NAME_FIELD]: [INV_DRUG_NAME, INV_DRUG_ID],
        [INV_DRUG_START_DAY_FIELD]: [INV_DRUG_START_DAY, INV_DRUG_START_DAY_CALCULATED],
        [INV_DRUG_END_DAY_FIELD]: [INV_DRUG_END_DAY, INV_DRUG_END_DAY_CALCULATED],
        [CON_MED_NAME_FIELD]: [CON_MED_TRT, CON_MED_DECOD],
        [CON_MED_START_DAY_FIELD]: [CON_MED_START_DAY, CON_MED_START_DAY_CALCULATED],
        [CON_MED_END_DAY_FIELD]: [CON_MED_END_DAY, CON_MED_END_DAY_CALCULATED],
        [MH_TERM_FIELD]: [MH_TERM, MH_DECOD_TERM]
    }

    globalStudyConfig = {
        [TRT_ARM_FIELD]: ACT_TRT_ARM,
        [AE_TERM_FIELD]: AE_TERM,
        [AE_START_DAY_FIELD]: AE_START_DAY,
        [AE_END_DAY_FIELD]: AE_END_DAY,
        [INV_DRUG_NAME_FIELD]: INV_DRUG_NAME,
        [INV_DRUG_START_DAY_FIELD]: INV_DRUG_START_DAY,
        [INV_DRUG_END_DAY_FIELD]: INV_DRUG_END_DAY,
        [CON_MED_NAME_FIELD]: CON_MED_TRT,
        [CON_MED_START_DAY_FIELD]: CON_MED_START_DAY,
        [CON_MED_END_DAY_FIELD]: CON_MED_END_DAY,
        [MH_TERM_FIELD]: MH_TERM
    }

    constructor(name) {
        super({});
        this.name = name;
    }

    createView(): void {
        let div = ui.divV([]);

        let inputs = ui.inputs([]);

        Object.keys(this.choices).forEach(field => {
            let fieldChoices = ui.choiceInput(`${field}`, this.choices[field][0], this.choices[field]);
            fieldChoices.onChanged((v) => {
                Object.keys(VIEWS_CONFIG).forEach(view => {
                    this.globalStudyConfig[field] = fieldChoices.value;
                    if (VIEWS_CONFIG[view][field]) {
                        VIEWS_CONFIG[view][field] = fieldChoices.value;
                        let obj = grok.shell.view(view) as any;
                        if (obj) {
                            if (obj.hasOwnProperty('loaded')) {
                                updateDivInnerHTML(obj.root, '');
                                obj.loaded = false;
                            } else {
                                const tableViewsParams = getTableViewsParams();
                                if(tableViewsParams[view]) {
                                    obj.close();
                                    const idx = VIEWS.findIndex(it => it.name === view);
                                    const tableView = createTableView(
                                        tableViewsParams[view].domainsAndColsToCheck,
                                        view,
                                        tableViewsParams[view].helpUrl,
                                        tableViewsParams[view].createViewHelper,
                                        tableViewsParams[view].paramsForHelper,
                                      );
                                    VIEWS[idx] = addView(tableView.view);
                                    VIEWS[idx].name = view;
                                }
                            }
                        }
                    }
                })
            })
            inputs.append(fieldChoices.root);
        });
        div.append(inputs);
        this.root.className = 'grok-view ui-box';
        this.root.append(div);
    }

}