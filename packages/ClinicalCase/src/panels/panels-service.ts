import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { CLIN_TRIAL_GOV_SEARCH, HttpService } from '../services/http.service';
import { addDataFromDmDomain, dictToString, getNullOrValue } from '../data-preparation/utils';
import { study } from '../clinical-study';
import { SEVERITY_COLOR_DICT, SUBJECT_ID, TREATMENT_ARM, AGE, SEX, RACE, CLINICAL_TRIAL_GOV_FIELDS } from '../constants';
import { AeBrowserView } from '../views/adverse-events-browser';
import { updateDivInnerHTML } from '../views/utils';
import $ from "cash-dom";
import { getSubjectDmData } from '../data-preparation/data-preparation';

export async function createPropertyPanel(viewClass: any) {
    switch (viewClass.name) {
        case 'Summary': {
            grok.shell.o = await summaryPanel(viewClass.studyId);
            break;
        }
        case 'Timelines': {
            grok.shell.o = await timelinesPanel(viewClass.resultTables, viewClass.selectedDataframes);
            break;
        }
        case 'AE Browser': {
            grok.shell.o = await aeBrowserPanel(viewClass);
        }
        default: {
            break;
        }
    }
}

export async function aeBrowserPanel(view: AeBrowserView) {
    if(view.aeToSelect.currentRowIdx !== -1){
        let subjId = view.aeToSelect.get(SUBJECT_ID, view.aeToSelect.currentRowIdx);
        let title = ui.tooltip.bind(ui.h1(subjId), dictToString(getSubjectDmData(subjId, [AGE, SEX, RACE, TREATMENT_ARM])));
        let description = ui.divH([view.aeToSelect.get('AETERM', view.aeToSelect.currentRowIdx)]);
        let severity = view.aeToSelect.get('AESEV', view.aeToSelect.currentRowIdx);
        let severityStyle = {style:{
            color: `${SEVERITY_COLOR_DICT[severity.toUpperCase()]}`,
              marginRight: '5px',
              fontWeight: 'bold'
        }}
        description.prepend(ui.divText(severity, severityStyle));
        let startEndDays = ui.tooltip.bind(
            ui.label(`${getNullOrValue(view.aeToSelect, 'AESTDY', view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, 'AEENDY', view.aeToSelect.currentRowIdx)}`), 
            `${getNullOrValue(view.aeToSelect, 'AESTDTC', view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, 'AEENDTC', view.aeToSelect.currentRowIdx)}`);
        let summary = ui.div([
          title,
          description,
          ui.inlineText(['Days ', startEndDays])
        ])
        let daysInput = ui.intInput('Days prior AE', view.daysPriorAe);
        daysInput.onChanged((v) => {
            view.daysPriorAe = daysInput.value;
            view.updateDomains();
            updateAccordion();
          });
        let accDiv = ui.div();
        let updateAccordion = () => {
            let acc = ui.accordion('ae-browser-panel');
            view.domains.concat(view.selectedAdditionalDomains).forEach(it => {
                const rowNum = view[ it ].rowCount === 1 && view[ it ].getCol(SUBJECT_ID).isNone(0) ? 0 : view[ it ].rowCount
                acc.addPane(`${it} (${rowNum})`, () => {
                    if (it) {
                        return ui.div(view[ it ].plot.grid().root)
                    }
                });
            })
            updateDivInnerHTML(accDiv, acc.root);
        }
        let addButton =  ui.icons.add(() => {
            let domainsMultiChoices = ui.multiChoiceInput('', view.selectedAdditionalDomains, view.additionalDomains)
            ui.dialog({ title: 'Select domains' })
              .add(ui.div([ domainsMultiChoices ]))
              .onOK(() => {
                view.selectedAdditionalDomains = domainsMultiChoices.value;
                view.updateDomains();
                updateAccordion();
              })
              //@ts-ignore
            .show({centerAt: addButton});
        });
        updateAccordion();
        return ui.div([
            summary,
            daysInput,
            accDiv,
            addButton
        ]);
    }
}

export async function summaryPanel(studyId: string) {
    const httpService = new HttpService();
    let clinTrialsGovInfo = await httpService.getStudyData('R01NS050536', Object.keys(CLINICAL_TRIAL_GOV_FIELDS));
    //let clinTrialsGovInfo = await httpService.getStudyData(studyId, Object.keys(CLINICAL_TRIAL_GOV_FIELDS));
    const summaryDict = {};
    Object.keys(clinTrialsGovInfo).forEach(key => {
        summaryDict[ CLINICAL_TRIAL_GOV_FIELDS[ key ] ] = clinTrialsGovInfo[ key ];
    })
    let studyLink = `${CLIN_TRIAL_GOV_SEARCH}${summaryDict[ 'NCT ID' ]}`
    summaryDict[ `Study link` ] = ui.link('Go to study page', () => { window.open(studyLink, '_blank').focus(); })
    let acc = ui.accordion('summary-panel');
    //acc.addPane(`${studyId}`, () => ui.tableFromMap(summaryDict), true);
    acc.addPane(`${studyId}`, () => {
        let info = ui.div();
        Object.keys(summaryDict).forEach(key => info.appendChild(
            ui.divH([
                ui.div(key, { style: { minWidth: '100px', fontWeight: 'bold' } }),
                ui.div(summaryDict[ key ])
            ])
        ))
        return info;
    });
    return acc.root;
}


export async function timelinesPanel(timelinesDf: DG.DataFrame, domains: string[]) {
    const selectedInd = timelinesDf.selection.getSelectedIndexes();
    if (!selectedInd.length) {
        if (domains.includes('ae')) {
            const aeWithArm = addDataFromDmDomain(timelinesDf.clone(), study.domains.dm, timelinesDf.columns.names(), [TREATMENT_ARM], 'key');
            const aeNumberByArm = aeWithArm.groupBy([ TREATMENT_ARM ]).where('domain = ae').count().aggregate();
            const subjNumberByArm = study.domains.dm.groupBy([ TREATMENT_ARM ]).uniqueCount(SUBJECT_ID).aggregate();
            const aeNumByArmDict = {};
            for (let i = 0; i < aeNumberByArm.rowCount; i++) {
                aeNumByArmDict[ aeNumberByArm.get(TREATMENT_ARM, i) ] = aeNumberByArm.get('count', i) /
                subjNumberByArm.
                groupBy([ TREATMENT_ARM, `unique(${SUBJECT_ID})` ])
                .where(`${TREATMENT_ARM} = ${aeNumberByArm.get(TREATMENT_ARM, i)}`)
                .aggregate()
                .get(`unique(${SUBJECT_ID})`, 0);
            }

            const aeTop5Dict = {};
            const aeTop5 = aeWithArm.groupBy([ TREATMENT_ARM, 'event' ]).where('domain = ae').count().aggregate();
            const order = aeTop5.getSortedOrder([ TREATMENT_ARM, 'event' ] as any);
            order.forEach(item => {
                const arm = aeTop5.get(TREATMENT_ARM, item);
                if (Object.keys(aeTop5Dict).includes(arm) && aeTop5Dict[ arm ].length < 5) {
                    aeTop5Dict[ arm ].push(aeTop5.get('event', item));
                } else {
                    aeTop5Dict[ arm ] = [ aeTop5.get('event', item) ]
                }
            })
            Object.keys(aeTop5Dict).forEach(key => {
                aeTop5Dict[ key ] = aeTop5Dict[ key ].join(',')
            })

            let root = ui.div();
            root.appendChild(ui.h2('Average AE number per subject'));
            root.appendChild(ui.tableFromMap(aeNumByArmDict));
            root.appendChild(ui.h2('Most frequent AEs'));
            root.appendChild(ui.tableFromMap(aeTop5Dict));
            return root;
        }
    } else {
        const eventArray = [];
        let summary;
        selectedInd.forEach((item, index) => {
            if (index === 0) {
                summary = { 'Subject ID': timelinesDf.get('key', item), 'Domain': timelinesDf.get('domain', item) }
            }
            const eventName = getNullOrValue(timelinesDf, 'event', item);
            const eventStart = getNullOrValue(timelinesDf, 'start', item);
            const eventEnd = getNullOrValue(timelinesDf, 'end', item);
            eventArray.push({
                event: eventName,
                start: eventStart,
                end: eventEnd
            })
        })
        let root = ui.div();
        root.appendChild(ui.tableFromMap(summary));
        const divsArray = [];
        eventArray.forEach(it => {
            divsArray.push(ui.div(ui.tableFromMap(it)));
        })
        root.appendChild(ui.divV(divsArray));
        return root;

    }

}