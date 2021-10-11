import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { CLIN_TRIAL_GOV_SEARCH, HttpService } from '../services/http.service';
import { addDataFromDmDomain, getNullOrValue } from '../data-preparation/utils';
import { study } from '../clinical-study';
import { SUBJECT_ID, TREATMENT_ARM } from '../constants';
import { AeBrowserView } from '../views/adverse-events-browser';

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
    let acc = ui.accordion();
    view.domains.forEach(it => {
        acc.addPane(it, () => {
            if (it) {
                return ui.div(view[ it ].plot.grid().root)
            }
        });
    })
    return acc.root;
}

export async function summaryPanel(studyId: string) {
    const httpService = new HttpService();
    let clinTrialsGovInfo = await httpService.getStudyData('R01NS050536');
    //let clinTrialsGovInfo = await httpService.getStudyData(studyId);
    let studyLink = `${CLIN_TRIAL_GOV_SEARCH}${clinTrialsGovInfo[ 'NCTId' ]}`
    clinTrialsGovInfo[ `Study link` ] = ui.link('Go to study page', () => { window.open(studyLink, '_blank').focus(); })
    let root = ui.div();
    root.appendChild(ui.h2(`Study ${studyId} summary`));
    root.appendChild(ui.tableFromMap(clinTrialsGovInfo));
    return root;
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