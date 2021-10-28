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
import { Accordion } from 'datagrok-api/dg';
import { AEBrowserHelper } from '../helpers/ae-browser-helper';

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

export async function aeBrowserPanel(view: AEBrowserHelper) {

    let panelDiv = ui.div();
    
    let accae = ui.accordion('ae-browser-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';

    if(view.aeToSelect.currentRowIdx !== -1){
        let subjId = view.aeToSelect.get(SUBJECT_ID, view.aeToSelect.currentRowIdx);
        let title = ui.tooltip.bind(ui.label(subjId), dictToString(getSubjectDmData(subjId, [AGE, SEX, RACE, TREATMENT_ARM])));
        let description = ui.divH([ui.divText(String(view.aeToSelect.get('AETERM', view.aeToSelect.currentRowIdx).toLowerCase()))]);
        let severity = view.aeToSelect.get('AESEV', view.aeToSelect.currentRowIdx);
        let severityStyle = {style:{
            color: `${SEVERITY_COLOR_DICT[severity.toUpperCase()]}`,
              marginRight: '5px',
              fontWeight: 'bold'
        }}

        accae.addTitle(ui.span([accIcon, title]));

        description.prepend(ui.divText(severity, severityStyle));
        let startEndDays = ui.tooltip.bind(
            ui.label(`${getNullOrValue(view.aeToSelect, 'AESTDY', view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, 'AEENDY', view.aeToSelect.currentRowIdx)}`), 
            `${getNullOrValue(view.aeToSelect, 'AESTDTC', view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, 'AEENDTC', view.aeToSelect.currentRowIdx)}`);
       

        let daysInput = ui.intInput('Prior AE', view.daysPriorAe);
        daysInput.onChanged((v) => {
            view.daysPriorAe = daysInput.value;
            view.updateDomains();
            updateAccordion();
          });
          startEndDays.innerHTML = 'Days ' + startEndDays.innerHTML;
          startEndDays.style.marginTop = '5px';
          //@ts-ignore
          accae.header = ui.div([
            description,
            startEndDays,
            //@ts-ignore
            ui.divH([ui.divText('Days prior AE'), daysInput.input], {style: {alignItems: 'center', gap: '5px'}}),
          ])

        let getPaneContent = (it, rowNum) => {
            if (it) {
                if(!rowNum){
                    return ui.divText('No records found');
                } else {
                    let grid = view[ it ].plot.grid();
                    if(rowNum < 7){
                        grid.root.style.maxHeight = rowNum < 4 ? '100px' : '150px' ;
                    }
                    grid.root.style.width = '250px';
                    return ui.div(grid.root);
                }
            }
        }

        let updateAccordion = () => {
            view.domains.concat(view.selectedAdditionalDomains).forEach(it => {
                const rowNum = view[ it ].rowCount === 1 && view[ it ].getCol(SUBJECT_ID).isNone(0) ? 0 : view[ it ].rowCount
                let pane = accae.getPane(`${it}`);
                //@ts-ignore
                updateDivInnerHTML(pane.root.lastChild, getPaneContent(it, rowNum));
                //@ts-ignore
                pane.root.firstChild.lastChild.innerText = rowNum;
            })
        }

        let createAccordion = () => {
            view.domains.concat(view.selectedAdditionalDomains).forEach(it => {
                const rowNum = view[ it ].rowCount === 1 && view[ it ].getCol(SUBJECT_ID).isNone(0) ? 0 : view[ it ].rowCount
                accae.addCountPane(`${it}`, () => getPaneContent(it, rowNum), () => rowNum);
                let panel = accae.getPane(`${it}`);
                //@ts-ignore
                $(panel.root).css('display', 'flex');
                //@ts-ignore
                $(panel.root).css('opacity', '1');
            })
        }
        
        let addButton =  ui.button(ui.icons.add(() => {}), ()=>{
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
        createAccordion();
        
        panelDiv.append(accae.root);
        panelDiv.append(addButton);

        return panelDiv;
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
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';

    acc.addTitle(ui.span([accIcon, ui.label(`${studyId}`)]));
    let acctable = ui.tableFromMap(summaryDict);
    acc.addPane('General', () => {
        $(acctable).find('tr').css('vertical-align', 'top');
        $(acctable).find('td').css('padding-bottom', '10px');
        $(acctable).find('.d4-entity-list>span').css('margin', '0px');
        return acctable
    }, true)

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
                aeTop5Dict[ key ] = aeTop5Dict[ key ].join(', ')
            })

            let acc = ui.accordion('timelines-panel')
            let accIcon = ui.element('i');
            accIcon.className = 'grok-icon svg-icon svg-view-layout';

            let avarageae = ui.tableFromMap(aeNumByArmDict);
            let mostae = ui.tableFromMap(aeTop5Dict);

            acc.addTitle(ui.span([accIcon, ui.label('Timelines')]));
            acc.addPane('Average AE', ()=>{
                $(avarageae).find('tr').css('vertical-align', 'top');
                $(avarageae).find('td').css('padding-bottom', '10px');
                $(avarageae).find('.d4-entity-list>span').css('margin', '0px');
                return avarageae
            }, true);

            acc.addPane('Frequents AEs', ()=>{
                $(mostae).find('tr').css('vertical-align', 'top');
                $(mostae).find('td').css('padding-bottom', '10px');
                $(mostae).find('.d4-entity-list>span').css('margin', '0px');
                return mostae
            });
            return acc.root;
        }
    } else {
        const eventArray = [];
        let acc2 = ui.accordion('timelines-patient-panel');
        let accIcon = ui.element('i');
        accIcon.className = 'grok-icon svg-icon svg-view-layout';

        let summary;
        selectedInd.forEach((item, index) => {
            if (index === 0) {
                summary = { 'Subject ID': timelinesDf.get('key', item), 'Domain': timelinesDf.get('domain', item) }
            }
            const eventName = String(getNullOrValue(timelinesDf, 'event', item)).toLowerCase();
            const eventStart = getNullOrValue(timelinesDf, 'start', item);
            const eventEnd = getNullOrValue(timelinesDf, 'end', item);
            eventArray.push({
                Event: eventName,
                Start: eventStart,
                End: eventEnd
            })
        })
        acc2.addTitle(ui.span([accIcon, ui.label(`${Object.values(summary)[0]}`)]));
        
        let summarytable = ui.tableFromMap(summary);
        acc2.addPane('General', ()=>{
            $(summarytable).find('tr').css('vertical-align', 'top');
            $(summarytable).find('td').css('padding-bottom', '10px');
            $(summarytable).find('.d4-entity-list>span').css('margin', '0px');
            return summarytable
        },true);
        
        const divsArray = [];
        eventArray.forEach(it => {
            let table = ui.tableFromMap(it);
            divsArray.push(ui.tableFromMap(it));
        })
        acc2.addPane('Events', ()=>{
            $(divsArray).find('.d4-entity-list>span').css('margin', '0px');
            return ui.divV(divsArray)
        });

        return acc2.root;

    }

}