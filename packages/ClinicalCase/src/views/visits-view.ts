import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains, updateDivInnerHTML } from './utils';
import { requiredColumnsByView } from '../constants';
import { createPivotedDataframe, filterNulls, getUniqueValues, getVisitNamesAndDays, addDataFromDmDomain } from '../data-preparation/utils';
import { AE_START_DAY, CON_MED_START_DAY, INV_DRUG_NAME, LAB_RES_N, SUBJECT_ID, SUBJ_REF_ENDT, TREATMENT_ARM, VISIT_DAY, VISIT_NAME, VISIT_START_DATE, VS_RES_N } from '../columns-constants';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';
import { createPropertyPanel } from '../panels/panels-service';
import { _package } from '../package';

export class VisitsView extends DG.ViewBase implements ILazyLoading {

    tv: DG.DataFrame;
    sv: DG.DataFrame;
    pivotedSv: DG.DataFrame;
    sortedVisitNamesAndDays: any;
    sortedVisitNames: any;
    patientVisit = new PatientVisit();
    studyVisit = new StudyVisit();
    totalVisits = {};
    proceduresAtVisit = { 'lb': {column: LAB_RES_N}, 'ex': {column: INV_DRUG_NAME}, 'vs': {column: VS_RES_N} };
    eventsSinceLastVisit = { 'ae': {column: AE_START_DAY}, 'cm': {column: CON_MED_START_DAY} };
    subjSet = new Set();
    existingDomains: string[];
    selectedDomain: string;
    selectedDomains: string[];
    heatMap: DG.Viewer;
    visitsGrid: DG.Grid;
    div = ui.box();
    ribbonDiv = ui.div();
    gridRibbonDiv: any;
    heatMapdomainChoices: any;
    selectedDomainsDiv = ui.div();

    constructor(name) {
        super({});
        this.name = name;
        this.helpUrl = `${_package.webRoot}/views_help/visits.md`;
    }

    loaded: boolean;

    load(): void {
        checkMissingDomains(requiredColumnsByView[this.name], this);
    }

    createView(): void {
        this.existingDomains = Object.keys(this.proceduresAtVisit)
            .concat(Object.keys(this.eventsSinceLastVisit))
            .filter(it => study.domains[it] !== null);
        this.assignColorsToDomains();
        this.tv = study.domains.tv.clone();
        this.sv = study.domains.sv.clone();
        filterNulls(this.sv, VISIT_DAY);
        this.pivotedSv = createPivotedDataframe(this.sv, [SUBJECT_ID], VISIT_NAME, VISIT_START_DATE, []);
        this.pivotedSv.columns.names().forEach(col => {
            if (this.pivotedSv.getCol(col).name !== SUBJECT_ID) {
                this.pivotedSv.getCol(col).tags.format = 'yyyy-MM-dd';
            }
            this.pivotedSv.getCol(col).name = col.replace(` avg(${VISIT_START_DATE})`, '');
        });
        this.sortedVisitNamesAndDays = getVisitNamesAndDays(this.tv, true);
        this.sortedVisitNames = this.sortedVisitNamesAndDays.map(it => it.name);
        let missingCols = this.sortedVisitNames.filter(it => !this.pivotedSv.columns.names().includes(it));
        missingCols.forEach(it => {
            this.pivotedSv.columns.addNewDateTime(it).init((i) => null);
        })
        this.pivotedSv = addDataFromDmDomain(this.pivotedSv, study.domains.dm, this.pivotedSv.columns.names(), [TREATMENT_ARM]);
        this.pivotedSv = this.pivotedSv.clone(null, [SUBJECT_ID, TREATMENT_ARM].concat(this.sortedVisitNames));
        this.pivotedSv.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
                this.createVisitPropertyPanel(this.pivotedSv);

            }, 100);

        });
        this.createTotalVisits();
        this.visitsGrid = this.pivotedSv.plot.grid();
        this.renderVisitCell();
        let switchGrid = this.createSwitchGridInput();
        this.createHeatMapRibbon();
        this.createGridRibbon();

        this.heatMap = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, this.createHeatMapDf(), {
            'colorCoding': 'All',
        });

        updateDivInnerHTML(this.div, this.visitsGrid.root);
        updateDivInnerHTML(this.ribbonDiv, this.gridRibbonDiv);

        this.root.className = 'grok-view ui-box';
        this.root.append(this.div);

        this.setRibbonPanels(
            [ [switchGrid.root], [this.ribbonDiv] ] ,
          );
    }

    private assignColorsToDomains() {
        this.existingDomains.forEach((domain, index) => {
            this.proceduresAtVisit[domain] ?
                this.proceduresAtVisit[domain]['color'] = DG.Color.toRgb(DG.Color.categoricalPalette[index]) :
                this.eventsSinceLastVisit[domain]['color'] = DG.Color.toRgb(DG.Color.categoricalPalette[index]);
        })
    }

    private createSwitchGridInput(){
        let switchGrid = ui.switchInput('Grid', true);
        switchGrid.onChanged((v) => {
            if(switchGrid.value){
                updateDivInnerHTML(this.div, this.visitsGrid.root);
                updateDivInnerHTML(this.ribbonDiv, this.gridRibbonDiv);
            } else {
                updateDivInnerHTML(this.div, this.heatMap.root);
                updateDivInnerHTML(this.ribbonDiv, this.heatMapdomainChoices.root);
            }
        });
        return switchGrid;
    }

    private createHeatMapRibbon() {
        this.selectedDomain = this.existingDomains[0];

        let domainChoices = ui.choiceInput('Domain', this.selectedDomain, this.existingDomains);
        domainChoices.onChanged((v) => {
            this.selectedDomain = domainChoices.value;
            this.heatMap.dataFrame = this.createHeatMapDf();
        });
        this.heatMapdomainChoices = domainChoices;
    }

    private createGridRibbon(){
        this.selectedDomains = this.existingDomains;
        let domainsButton =  ui.button('Domains', ()=>{
            //@ts-ignore
            let domainsMultiChoices = ui.multiChoiceInput('', this.selectedDomains, this.existingDomains);
            ui.dialog({ title: 'Select domains' })
              .add(ui.div([ domainsMultiChoices ]))
              .onOK(() => {
                this.selectedDomains = domainsMultiChoices.value;
                updateDivInnerHTML(this.selectedDomainsDiv, this.createSelectedDomainsDiv());
                this.visitsGrid.invalidate();
              })
              //@ts-ignore
            .show();
        });

        this.gridRibbonDiv = ui.divH([
            domainsButton,
            this.selectedDomainsDiv
        ]);
        updateDivInnerHTML(this.selectedDomainsDiv, this.createSelectedDomainsDiv());
    }


    private createSelectedDomainsDiv(){
        let selectedDomainsDiv = ui.divH([]);
        this.selectedDomains.forEach((it) => {
            let domainName = ui.divText(`${it} `);
            domainName.style.color = this.proceduresAtVisit[it] ? this.proceduresAtVisit[it].color : this.eventsSinceLastVisit[it].color;
            domainName.style.paddingRight = '7px';
            domainName.style.paddingTop = '10px';
            selectedDomainsDiv.append(domainName);
        });
       return selectedDomainsDiv;
    }


    private createHeatMapDf() {
        let df = DG.DataFrame.create();
        df.columns.addNewString(SUBJECT_ID);
        this.pivotedSv.columns.names().forEach(col => {
            if (col !== SUBJECT_ID && col !== TREATMENT_ARM) {
                df.columns.addNewInt(col);
            }
        });
        Array.from(this.subjSet).forEach(id => {
            df.rows.addNew();
            df.set(SUBJECT_ID, df.rowCount - 1, id);
            Object.keys(this.totalVisits).forEach(key => {
                df.set(key, df.rowCount - 1, this.totalVisits[key][id].eventsCount[this.selectedDomain]);
            })
        });
        df = addDataFromDmDomain(df, study.domains.dm, df.columns.names().filter(it => it !== TREATMENT_ARM), [TREATMENT_ARM]);
        df = df.clone(null, [SUBJECT_ID, TREATMENT_ARM].concat(this.sortedVisitNames));
        this.setColorPaletteForHeatMap(df);
        df.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
                this.createVisitPropertyPanel(df);

            }, 100);

        });
        return df;
    }

    private setColorPaletteForHeatMap(df: DG.DataFrame){
        df.columns.names().forEach(col => {
            if (col !== SUBJECT_ID && col !== TREATMENT_ARM) {
                df.col(col).tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
                df.col(col).tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.white}, ${DG.Color.blue}]`;
            }
        });
    }

    private createVisitPropertyPanel(df: DG.DataFrame){
        if (df.currentCol.name !== SUBJECT_ID && df.currentCol.name !== TREATMENT_ARM) {
            if (df.currentRowIdx === -1) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(df.currentCol.name);
                this.studyVisit.createStudyVisitPanel(currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
            } else {
                let subjId = df.get(SUBJECT_ID, df.currentRowIdx);
                let currentPatientVisit = this.totalVisits[df.currentCol.name][subjId]
                currentPatientVisit.updateSubjectVisit();
                createPropertyPanel(currentPatientVisit);
            }
        }

    }

    getCurrentAndPreviousVisits(colName: string) {
        let visit = this.sortedVisitNamesAndDays.filter(it => it.name === colName)[0];
        let visitIdx = this.sortedVisitNamesAndDays.findIndex(it => it.name === colName);
        let previousVisit = visitIdx > 0 ? this.sortedVisitNamesAndDays[visitIdx - 1] : null;
        return { current: visit, previous: previousVisit };
    }

    createTotalVisits() {
        let countDfs = this.datasetsWithNumberProceduresAtVisit();
        this.createInitialTotalVisits();
        this.updateProceduresAtVisitCount(countDfs);
        this.updateEventsSinceLastVisitCount();
    }

    private datasetsWithNumberProceduresAtVisit() {
        let countDfs = {};
        Object.keys(this.proceduresAtVisit).forEach(domain => {
            if (study.domains[domain]) {
                countDfs[domain] = study.domains[domain]
                    .groupBy([SUBJECT_ID, VISIT_NAME])
                    .count(this.proceduresAtVisit[domain].column)
                    .aggregate();
            }
        });
        return countDfs;
    }

    private createInitialTotalVisits() {
        this.pivotedSv.columns.names().forEach(colName => {
            if (colName !== SUBJECT_ID && colName !== TREATMENT_ARM) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(colName);
                this.totalVisits[colName] = {};

                this.pivotedSv.getCol(SUBJECT_ID).categories.forEach(subjId => {
                    this.totalVisits[colName][subjId] = new PatientVisit();
                    this.totalVisits[colName][subjId].currentSubjId = subjId;
                    this.totalVisits[colName][subjId].currentVisitDay = currentVisit.day;
                    this.totalVisits[colName][subjId].currentVisitName = currentVisit.name;
                    this.totalVisits[colName][subjId].previsousVisitDay = previousVisit ? previousVisit.day : null;
                });
            }

        });
    }

    private updateProceduresAtVisitCount(countDfs: any) {
        Object.keys(countDfs).forEach(domain => {
            for (let i = 0; i < countDfs[domain].rowCount; i++) {
                let visitName = countDfs[domain].get(VISIT_NAME, i);
                if (!this.sortedVisitNames.includes(visitName)) {
                    continue;
                }
                let subjId = countDfs[domain].get(SUBJECT_ID, i);
                this.subjSet.add(subjId);
                this.totalVisits[visitName][subjId].eventsCount[domain] = countDfs[domain].get(this.proceduresAtVisit[domain].column, i);
            }
        });
    }

    private updateEventsSinceLastVisitCount() {
        Object.keys(this.eventsSinceLastVisit).forEach(domain => {
            if (study.domains[domain]) {
                for (let i = 0; i < study.domains[domain].rowCount; i++) {
                    let startDay = study.domains[domain].get(this.eventsSinceLastVisit[domain].column, i);
                    let subjId = study.domains[domain].get(SUBJECT_ID, i);
                    this.subjSet.add(subjId);
                    for (let z = 0; z < this.sortedVisitNamesAndDays.length - 1; z++) {
                        if (startDay > this.sortedVisitNamesAndDays[z].day && startDay < this.sortedVisitNamesAndDays[z + 1].day) {
                            let visitName = this.sortedVisitNamesAndDays[z + 1].name;
                            if (this.totalVisits[visitName][subjId].eventsCount[domain]) {
                                this.totalVisits[visitName][subjId].eventsCount[domain] += 1;
                            } else {
                                this.totalVisits[visitName][subjId].eventsCount[domain] = 1;
                            }
                            break;
                        }
                    }
                }
            }
        })
    }

    private renderVisitCell() {
        this.visitsGrid.onCellRender.subscribe((args) => {
            let gc = args.cell;
            if (gc.isTableCell && gc.gridColumn.name !== SUBJECT_ID && gc.gridColumn.name !== TREATMENT_ARM) {
                let patientVisit = this.totalVisits[gc.gridColumn.name][this.visitsGrid.dataFrame.get(SUBJECT_ID, gc.tableRowIndex)];
                let delta = 30;
                let x = args.bounds.x + 10;
                let y = args.bounds.y + 15;
                this.selectedDomains.forEach((item) => {
                    if (patientVisit.eventsCount[item]) {
                        args.g.fillStyle = this.proceduresAtVisit[item] ? this.proceduresAtVisit[item].color : this.eventsSinceLastVisit[item].color;
                        args.g.fillText(patientVisit.eventsCount[item], x, y);
                    } else {
                        args.g.fillStyle = "rgba(128, 128, 128, 0.1)";
                        args.g.fillText('0', x, y);
                    }
                    x += delta;
                })
                args.preventDefault();
            }
        })

    }

}
