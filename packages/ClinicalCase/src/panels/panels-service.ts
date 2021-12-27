import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";
import { CLIN_TRIAL_GOV_SEARCH, HttpService } from '../services/http.service';
import { addDataFromDmDomain, dictToString, getNullOrValue } from '../data-preparation/utils';
import { study } from '../clinical-study';
import { SEVERITY_COLOR_DICT, CLINICAL_TRIAL_GOV_FIELDS } from '../constants';
import { SUBJECT_ID, TREATMENT_ARM, AGE, SEX, RACE, AE_TERM, AE_SEVERITY, AE_START_DAY, AE_END_DAY, AE_START_DATE, AE_END_DATE, DOMAIN, INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, CON_MED_DOSE, CON_MED_DOSE_UNITS, CON_MED_DOSE_FREQ, CON_MED_ROUTE, INV_DRUG_DOSE_FORM, INV_DRUG_DOSE_FREQ, INV_DRUG_ROUTE, AE_DECOD_TERM, CON_MED_NAME, LAB_TEST, LAB_RES_N, VS_TEST, VS_RES_N } from '../columns-constants';
import { updateDivInnerHTML } from '../views/utils';
import $ from "cash-dom";
import { getSubjectDmData } from '../data-preparation/data-preparation';
import { AEBrowserHelper } from '../helpers/ae-browser-helper';
import { LaboratoryView } from '../views/laboratory-view';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';

export async function createPropertyPanel(viewClass: any) {
    switch (viewClass.name) {
        case 'Summary': {
            grok.shell.o = await summaryPanel(viewClass.studyId);
            break;
        }
        case 'Timelines': {
            const panel = await timelinesPanel(viewClass.resultTables, viewClass.selectedDataframes, viewClass.aeBrowserHelper);
            if (panel) {
                grok.shell.o = panel;
            }
            break;
        }
        case 'AE Browser': {
            grok.shell.o = await aeBrowserPanel(viewClass);
        }
        case 'Laboratory': {
            grok.shell.o = await laboratoryPanel(viewClass);
            break;
        }
        case 'Patient Visit': {
            grok.shell.o = await patientVisitPanel(viewClass);
            break;
        }
        case 'Study Visit': {
            grok.shell.o = await studyVisitPanel(viewClass);
            break;
        }
        default: {
            break;
        }
    }
}

export async function studyVisitPanel(studyVisit: StudyVisit) {
    let panelDiv = ui.div();

    let acc = ui.accordion('patient-visit-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(studyVisit.currentVisitName)]));

    acc.addPane(`Summary`, () => {
        return ui.tableFromMap({
            'Visit day': studyVisit.currentVisitDay,
            'Total patients': studyVisit.totalPatients,
            'Min visit date': studyVisit.minVisitDate,
            'Max visit dae': studyVisit.maxVisitDate,
        })
    });

    let getRowNumber = (df) => {
        return df ? df.rowCount === 1 && df.getCol(SUBJECT_ID).isNone(0) ? 0 : df.rowCount : 0;
    }

    acc.addPane('Drug exposure', () => {
        if (!getRowNumber(studyVisit.exAtVisit)) {
            return ui.divText('No records found');
        }
        return DG.Viewer.fromType(DG.VIEWER.PIE_CHART, studyVisit.exAtVisit, {
            category: studyVisit.extrtWithDoseColName,
        }).root;
    })

    let createPane = (name, rowNum, df, splitCol) => {
        acc.addCountPane(`${name}`, () => getPaneContent(df, splitCol, rowNum), () => rowNum);
        let panel = acc.getPane(`${name}`);
        //@ts-ignore
        $(panel.root).css('display', 'flex');
        //@ts-ignore
        $(panel.root).css('opacity', '1');
    }

    let getPaneContent = (df, splitCol, rowNum) => {
            if (!rowNum) {
                return ui.divText('No records found');
            } else {
                return df.plot.bar({
                    split: splitCol,
                    style: 'dashboard',
                    legendVisibility: 'Never'
                }).root;
            }
    }

    let aeRowNum = getRowNumber(studyVisit.aeSincePreviusVisit);
    createPane('AEs since previous visit', aeRowNum, studyVisit.aeSincePreviusVisit, AE_DECOD_TERM);

    let cmRowNum = getRowNumber(studyVisit.conmedSincePreviusVisit);
    createPane('CONMEDs since previous visit', cmRowNum, studyVisit.conmedSincePreviusVisit, CON_MED_NAME);

    let createDistributionPane = (name, df, catCol, valCol) => {
        acc.addPane(name, () => {
            if (!getRowNumber(df)) {
                return ui.divText('No records found');
            }
            let categoriesAcc = ui.accordion();
            df.getCol(catCol).categories.forEach(cat => {
                categoriesAcc.addPane(`${cat}`, () => {
                    let valueDf = df.groupBy(df.columns.names())
                    .where(`${catCol} = ${cat}`)
                    .aggregate();
                    const plot = DG.Viewer.boxPlot(valueDf, {
                        value: `${valCol}`,
                        category: `${catCol}`,
                        labelOrientation: 'Horz',
                        showCategorySelector: false,
                        showValueSelector: false,
                        showPValue: true
                      });
                    return plot.root;
                })
            })
            return categoriesAcc.root;   
        })
    } 

    createDistributionPane('Laboratory', studyVisit.lbAtVisit, LAB_TEST, LAB_RES_N);
    createDistributionPane('Vital signs', studyVisit.vsAtVisit, VS_TEST, VS_RES_N);

    panelDiv.append(acc.root);

    return panelDiv;

}

export async function patientVisitPanel(patientVisit: PatientVisit) {
    let panelDiv = ui.div();

    let acc = ui.accordion('patient-visit-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acc.addTitle(ui.span([accIcon, ui.label(patientVisit.currentSubjId)]));

    let getPaneContent = (it, rowNum) => {
        if (it) {
            if (!rowNum) {
                return ui.divText('No records found');
            } else {
                let grid = patientVisit[it].plot.grid();
                if (rowNum < 7) {
                    grid.root.style.maxHeight = rowNum < 4 ? '100px' : '150px';
                }
                grid.root.style.width = '250px';
                return ui.div(grid.root);
            }
        }
    }

    let createPane = (it, name, rowNum) => {
        acc.addCountPane(`${name}`, () => getPaneContent(it, rowNum), () => rowNum);
        let panel = acc.getPane(`${name}`);
        //@ts-ignore
        $(panel.root).css('display', 'flex');
        //@ts-ignore
        $(panel.root).css('opacity', '1');
    }


    let createAccordion = () => {
        Object.keys(patientVisit.domainsNamesDict).forEach(key => {
            let domain = patientVisit.domainsNamesDict[key];
            const rowNum = patientVisit[domain] ?
                patientVisit[domain].rowCount === 1 && patientVisit[domain].getCol(SUBJECT_ID).isNone(0) ?
                    0 : patientVisit[domain].rowCount : 0
            createPane(domain, key, rowNum);
        })
    }

    createAccordion();

    panelDiv.append(acc.root);

    return panelDiv;

}

export async function laboratoryPanel(view: LaboratoryView){
    let panelDiv = ui.div();  
    let acclab = ui.accordion('laboratory-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';
    acclab.addTitle(ui.span([accIcon, 'Laboratory']));

    let altChoices = ui.choiceInput('ALT', view.selectedALT, view.uniqueLabValues);
    altChoices.onChanged((v) => {
      view.selectedALT = altChoices.value;
      view.updateHysLawScatterPlot();
    });
    //@ts-ignore
    altChoices.input.style.width = '150px';

    let astChoices = ui.choiceInput('AST', view.selectedAST, view.uniqueLabValues);
    astChoices.onChanged((v) => {
      view.selectedAST = astChoices.value;
      view.updateHysLawScatterPlot();
    });
    //@ts-ignore
    astChoices.input.style.width = '150px';

    let blnChoices = ui.choiceInput('BLN', view.selectedBLN, view.uniqueLabValues);
    blnChoices.onChanged((v) => {
        view.selectedBLN = blnChoices.value;
        view.updateHysLawScatterPlot();
    });
    //@ts-ignore
    blnChoices.input.style.width = '150px';
    acclab.addPane('Hy\'s law', () => ui.divV([ altChoices.root, astChoices.root, blnChoices.root ]), true);

    let labValueChoices = ui.choiceInput('Value', view.selectedLabBlEp, view.uniqueLabValues);
    let blVisitChoices = ui.choiceInput('BL', view.selectedBl, view.uniqueVisits);
    let epVisitChoices = ui.choiceInput('EP', view.selectedEp, view.uniqueVisits);
    labValueChoices.onChanged((v) => {
        view.selectedLabBlEp = labValueChoices.value;
        view.updateBaselineEndpointPlot();
    });
    blVisitChoices.onChanged((v) => {
        view.selectedBl = blVisitChoices.value;
        view.updateBaselineEndpointPlot();
    });
    epVisitChoices.onChanged((v) => {
        view.selectedEp = epVisitChoices.value;
        view.updateBaselineEndpointPlot();
    });
    //@ts-ignore
    labValueChoices.input.style.width = '150px';
    acclab.addPane('Baseline endpoint', () => ui.divV([ labValueChoices.root, blVisitChoices.root, epVisitChoices.root ]), true);


    let distrChoices = ui.choiceInput('Value', view.selectedLabDistr, view.uniqueLabValues);
    let treatmentArmsChoices = ui.choiceInput('Treatment arm', view.selectedArm, view.uniqueTreatmentArms);
    distrChoices.onChanged((v) => {
      view.selectedLabDistr = distrChoices.value;
      view.updateDistributionPlot();
    });
    treatmentArmsChoices.onChanged((v) => {
      view.selectedArm = treatmentArmsChoices.value;
      view.updateDistributionPlot();
    });
    //@ts-ignore
    distrChoices.input.style.width = '150px';
    acclab.addPane('Laboratory distribution', () => ui.divV([ distrChoices.root, treatmentArmsChoices.root ]), true);
    panelDiv.append(acclab.root);
    return panelDiv;
}

export async function aeBrowserPanel(view: AEBrowserHelper) {

    let panelDiv = ui.div();
    
    let accae = ui.accordion('ae-browser-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';

    if(view.aeToSelect.currentRowIdx !== -1){
        let subjId = view.aeToSelect.get(SUBJECT_ID, view.aeToSelect.currentRowIdx);
        let title = ui.tooltip.bind(ui.label(subjId), dictToString(getSubjectDmData(subjId, [AGE, SEX, RACE, TREATMENT_ARM])));
        let description = ui.divH([ui.divText(String(view.aeToSelect.get(AE_TERM, view.aeToSelect.currentRowIdx).toLowerCase()))]);
        let severity = view.aeToSelect.get(AE_SEVERITY, view.aeToSelect.currentRowIdx);
        let severityStyle = {style:{
            color: `${SEVERITY_COLOR_DICT[severity.toUpperCase()]}`,
              marginRight: '5px',
              fontWeight: 'bold'
        }}

        accae.addTitle(ui.span([accIcon, title]));

        description.prepend(ui.divText(severity, severityStyle));
        let startEndDays = ui.tooltip.bind(
            ui.label(`${getNullOrValue(view.aeToSelect, AE_START_DAY, view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, AE_END_DAY, view.aeToSelect.currentRowIdx)}`), 
            `${getNullOrValue(view.aeToSelect, AE_START_DATE, view.aeToSelect.currentRowIdx)} - ${getNullOrValue(view.aeToSelect, AE_END_DATE, view.aeToSelect.currentRowIdx)}`);
       

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

        let createPane = (it, rowNum) => {
            accae.addCountPane(`${it}`, () => getPaneContent(it, rowNum), () => rowNum);
            let panel = accae.getPane(`${it}`);
            //@ts-ignore
            $(panel.root).css('display', 'flex');
            //@ts-ignore
            $(panel.root).css('opacity', '1');
        }

        let updateAccordion = () => {
            const totalDomains =  view.domains.concat(view.selectedAdditionalDomains);
            const panesToRemove = accae.panes.filter(it => !totalDomains.includes(it.name));
            panesToRemove.forEach(it => accae.removePane(it));
            totalDomains.forEach(it => {
                const rowNum = view[it].rowCount === 1 && view[it].getCol(SUBJECT_ID).isNone(0) ? 0 : view[it].rowCount
                let pane = accae.getPane(`${it}`);
                if (pane) {
                    //@ts-ignore
                    updateDivInnerHTML(pane.root.lastChild, getPaneContent(it, rowNum));
                    //@ts-ignore
                    pane.root.firstChild.lastChild.innerText = rowNum;
                } else {
                    createPane(it, rowNum);
                }
            })
        }

        let createAccordion = () => {
            view.domains.concat(view.selectedAdditionalDomains).forEach(it => {
                const rowNum = view[ it ].rowCount === 1 && view[ it ].getCol(SUBJECT_ID).isNone(0) ? 0 : view[ it ].rowCount
                createPane(it, rowNum);
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
    //let clinTrialsGovInfo = await httpService.getStudyData('R01NS050536', Object.keys(CLINICAL_TRIAL_GOV_FIELDS));
    let clinTrialsGovInfo = await httpService.getStudyData(studyId, Object.keys(CLINICAL_TRIAL_GOV_FIELDS));

    let acc = ui.accordion('summary-panel');
    let accIcon = ui.element('i');
    accIcon.className = 'grok-icon svg-icon svg-view-layout';

    acc.addTitle(ui.span([accIcon, ui.label(`${studyId}`)]));

    if (clinTrialsGovInfo) {
        const summaryDict = {};
        Object.keys(clinTrialsGovInfo).forEach(key => {
            summaryDict[CLINICAL_TRIAL_GOV_FIELDS[key]] = clinTrialsGovInfo[key];
        })
        let studyLink = `${CLIN_TRIAL_GOV_SEARCH}${summaryDict['NCT ID']}`
        summaryDict[`Study link`] = ui.link('Go to study page', () => { window.open(studyLink, '_blank').focus(); })

        let acctable = ui.tableFromMap(summaryDict);
        acc.addPane('General', () => {
            $(acctable).find('tr').css('vertical-align', 'top');
            $(acctable).find('td').css('padding-bottom', '10px');
            $(acctable).find('.d4-entity-list>span').css('margin', '0px');
            return acctable
        }, true)
    } else {
        acc.addPane('General', () => {
            return ui.divText('Study not found on clinicaltrials.gov')
        }, true)
    }

    return acc.root;
}


export async function timelinesPanel(timelinesDf: DG.DataFrame, domains: string[], aeBrowserHelper: AEBrowserHelper) {
    const selectedInd = timelinesDf.selection.getSelectedIndexes();

    if (!selectedInd.length) {
        if (domains.includes('ae')) {
            const aeWithArm = addDataFromDmDomain(timelinesDf.clone(), study.domains.dm, timelinesDf.columns.names(), [TREATMENT_ARM], 'key');
            const aeNumberByArm = aeWithArm.groupBy([TREATMENT_ARM]).where(`${DOMAIN} = ae`).count().aggregate();
            const subjNumberByArm = study.domains.dm.groupBy([TREATMENT_ARM]).uniqueCount(SUBJECT_ID).aggregate();
            const aeNumByArmDict = {};
            for (let i = 0; i < aeNumberByArm.rowCount; i++) {
                aeNumByArmDict[aeNumberByArm.get(TREATMENT_ARM, i)] = aeNumberByArm.get('count', i) /
                    subjNumberByArm.
                        groupBy([TREATMENT_ARM, `unique(${SUBJECT_ID})`])
                        .where(`${TREATMENT_ARM} = ${aeNumberByArm.get(TREATMENT_ARM, i)}`)
                        .aggregate()
                        .get(`unique(${SUBJECT_ID})`, 0);
            }

            const aeTop5Dict = {};
            const aeTop5 = aeWithArm.groupBy([TREATMENT_ARM, 'event']).where(`${DOMAIN} = ae`).count().aggregate();
            const order = aeTop5.getSortedOrder([TREATMENT_ARM, 'event'] as any);
            order.forEach(item => {
                const arm = aeTop5.get(TREATMENT_ARM, item);
                if (Object.keys(aeTop5Dict).includes(arm) && aeTop5Dict[arm].length < 5) {
                    aeTop5Dict[arm].push(aeTop5.get('event', item));
                } else {
                    aeTop5Dict[arm] = [aeTop5.get('event', item)]
                }
            })

            Object.keys(aeTop5Dict).forEach(key => {
                aeTop5Dict[key] = aeTop5Dict[key].join(', ')
            })

            let acc = ui.accordion('timelines-panel')
            let accIcon = ui.element('i');
            accIcon.className = 'grok-icon svg-icon svg-view-layout';

            let avarageae = ui.tableFromMap(aeNumByArmDict);
            let mostae = ui.tableFromMap(aeTop5Dict);

            acc.addTitle(ui.span([accIcon, ui.label('Timelines')]));
            acc.addPane('Average AE per patient', () => {
                $(avarageae).find('tr').css('vertical-align', 'top');
                $(avarageae).find('td').css('padding-bottom', '10px');
                $(avarageae).find('.d4-entity-list>span').css('margin', '0px');
                return avarageae
            }, true);

            acc.addPane('Frequents AEs', () => {
                $(mostae).find('tr').css('vertical-align', 'top');
                $(mostae).find('td').css('padding-bottom', '10px');
                $(mostae).find('.d4-entity-list>span').css('margin', '0px');
                return mostae
            });
            return acc.root;
        }
    } else {
        const eventArray = [];
        const eventIndexesArray = [];
        let acc2 = ui.accordion('timelines-patient-panel');
        let accIcon = ui.element('i');
        accIcon.className = 'grok-icon svg-icon svg-view-layout';

        let domainAdditionalFields = {
            ae: { fields: [AE_SEVERITY], pos: 'before' },
            cm: { fields: [CON_MED_DOSE, CON_MED_DOSE_UNITS, CON_MED_DOSE_FREQ, CON_MED_ROUTE], pos: 'after' },
            ex: { fields: [INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, INV_DRUG_DOSE_FORM, INV_DRUG_DOSE_FREQ, INV_DRUG_ROUTE], pos: 'after' },
        }

        let switchToAEBrowserPanel = (aeRowNum) => {
            if (aeBrowserHelper.aeToSelect.currentRowIdx === aeRowNum) {
                aeBrowserHelper.createAEBrowserPanel();
            } else {
                //@ts-ignore
                aeBrowserHelper.aeToSelect.currentRow = aeRowNum;
            }
        }

        if (selectedInd.length === 1 && timelinesDf.get('domain', selectedInd[0]) === 'ae') {
            const aeRowNum = timelinesDf.get('rowNum', selectedInd[0]);
            switchToAEBrowserPanel(aeRowNum);
            return null;
        }

        selectedInd.forEach((item) => {
            const domain = timelinesDf.get('domain', item);
            const addDomainInfo = domainAdditionalFields[domain];
            let addInfoString = '';
            const index = timelinesDf.get('rowNum', item);
            addDomainInfo.fields.forEach(it => {
                if (study.domains[domain].columns.names().includes(it)) {
                    addInfoString += `${study.domains[domain].get(it, index)} `;
                }
            });
            const eventName = String(getNullOrValue(timelinesDf, 'event', item)).toLowerCase();
            const fullEventName = addDomainInfo.pos === 'before' ? `${addInfoString}${eventName}` : `${eventName} ${addInfoString}`;
            const eventStart = getNullOrValue(timelinesDf, 'start', item);
            const eventEnd = getNullOrValue(timelinesDf, 'end', item);
            eventArray.push({
                Domain: domain,
                Event: fullEventName,
                Days: `${eventStart} - ${eventEnd}`
            })
            eventIndexesArray.push(index);
        })
        acc2.addTitle(ui.span([accIcon, ui.label(`${timelinesDf.get('key', selectedInd[0])}`)]));

        const eventTable = DG.DataFrame.fromObjects(eventArray);
        const eventGrid = eventTable.plot.grid();
        eventGrid.columns.byName('domain').width = 55;
        let col = eventGrid.columns.byName('event');
        col.width = 170;
        col.cellType = 'html';

        eventGrid.onCellPrepare(function (gc) {
            if (gc.isTableCell && eventTable.get('Domain', gc.gridRow) === 'ae' && gc.gridColumn.name === 'Event') {
                let eventElement = ui.link(gc.cell.value, {}, '', { id: `${eventIndexesArray[gc.gridRow]}` });
                eventElement.addEventListener('click', (event) => {
                    switchToAEBrowserPanel(parseInt(eventElement.id));
                    event.stopPropagation();
                });
                gc.style.element = ui.div(eventElement, {style: {'white-space': 'nowrap'}});
            } else {
                gc.style.element = ui.divText(gc.cell.value, {style: {'white-space': 'nowrap'}});
            }
            gc.style.element.style.paddingTop = '7px';
            gc.style.element.style.paddingLeft = '7px';
            ui.tooltip.bind(gc.style.element, gc.cell.value);
        });

        acc2.addPane('Events', () => {
            return ui.div(eventGrid.root);
        });
        const accPane = acc2.getPane('Events').root.getElementsByClassName('d4-accordion-pane-content')[0] as HTMLElement;
        accPane.style.margin = '0px';
        accPane.style.paddingLeft = '0px';

        return acc2.root;
    }

}