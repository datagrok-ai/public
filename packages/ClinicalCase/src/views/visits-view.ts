import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { updateDivInnerHTML } from './utils';
import { createPivotedDataframe, filterNulls, getVisitNamesAndDays, addDataFromDmDomain, dataframeContentToRow } from '../data-preparation/utils';
import { AE_DECOD_TERM, AE_START_DAY, CON_MED_NAME, CON_MED_START_DAY, INV_DRUG_NAME, LAB_RES_N, LAB_TEST, SUBJECT_ID, TREATMENT_ARM, VISIT_DAY, VISIT_NAME, VISIT_START_DATE, VS_RES_N, VS_TEST } from '../columns-constants';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';
import { _package } from '../package';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';
import { DockManager } from 'datagrok-api/dg';

export class VisitsView extends ClinicalCaseViewBase {

    tv: DG.DataFrame;
    sv: DG.DataFrame;
    pivotedSv: DG.DataFrame;
    sortedVisitNamesAndDays: any;
    sortedVisitNames: any;
    patientVisit = new PatientVisit(study.domains);
    studyVisit = new StudyVisit(study.domains);
    totalVisits = {};
    proceduresAtVisit = { 'lb': { column: LAB_TEST}, 'ex': { column: INV_DRUG_NAME }, 'vs': { column: VS_TEST } };
    eventsSinceLastVisit = { 'ae': { column: AE_START_DAY }, 'cm': { column: CON_MED_START_DAY } };
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

    createView(): void {
        this.existingDomains = Object.keys(this.proceduresAtVisit)
            .concat(Object.keys(this.eventsSinceLastVisit))
            .filter(it => study.domains[it] !== null);
        this.assignColorsToDomains();
        this.tv = study.domains.tv.clone();
        this.sv = study.domains.sv.clone();
        filterNulls(this.sv, VISIT_DAY); //remove unshc visits
        this.pivotedSv = createPivotedDataframe(this.sv, [SUBJECT_ID], VISIT_NAME, VISIT_START_DATE, []);
        this.pivotedSv.columns.names().forEach(col => {
            if (this.pivotedSv.getCol(col).name !== VISIT_START_DATE) {
                this.pivotedSv.getCol(col).tags.format = 'yyyy-MM-dd';
            }
            this.pivotedSv.getCol(col).name = col.replace(` avg(${SUBJECT_ID})`, '');
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
            [[switchGrid.root], [this.ribbonDiv]] ,
        );
    }

    private assignColorsToDomains() {
        this.existingDomains.forEach((domain, index) => {
            this.proceduresAtVisit[domain] ?
                this.proceduresAtVisit[domain]['color'] = DG.Color.toRgb(DG.Color.categoricalPalette[index]) :
                this.eventsSinceLastVisit[domain]['color'] = DG.Color.toRgb(DG.Color.categoricalPalette[index]);
        })
    }

    private createSwitchGridInput() {
        let switchGrid = ui.switchInput('Grid', true);
        switchGrid.onChanged((v) => {
            if (switchGrid.value) {
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

    private createGridRibbon() {
        this.selectedDomains = this.existingDomains;
        let domainsButton = ui.button('Domains', () => {
            //@ts-ignore
            let domainsMultiChoices = ui.multiChoiceInput('', this.selectedDomains, this.existingDomains);
            ui.dialog({ title: 'Select domains' })
                .add(ui.div([domainsMultiChoices]))
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


    private createSelectedDomainsDiv() {
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
        typeof(study.domains.dm.get(SUBJECT_ID, 0)) === 'number' ? df.columns.addNewInt(SUBJECT_ID) : df.columns.addNewString(SUBJECT_ID);
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
        console.log(df.columns.names())
        df = df.clone(null, [SUBJECT_ID, TREATMENT_ARM].concat(this.sortedVisitNames));
        this.setColorPaletteForHeatMap(df);
        df.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
                this.createVisitPropertyPanel(df);

            }, 100);

        });
        return df;
    }

    private setColorPaletteForHeatMap(df: DG.DataFrame) {
        df.columns.names().forEach(col => {
            if (col !== SUBJECT_ID && col !== TREATMENT_ARM) {
                df.col(col).tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
                df.col(col).tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.white}, ${DG.Color.blue}]`;
            }
        });
    }

    private async createVisitPropertyPanel(df: DG.DataFrame) {
        if (df.currentCol.name !== SUBJECT_ID && df.currentCol.name !== TREATMENT_ARM) {
            if (df.currentRowIdx === -1) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(df.currentCol.name);
                this.studyVisit.updateStudyVisit(currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
                grok.shell.o = await this.studyVisitPanel();
            } else {
                let subjId = df.get(SUBJECT_ID, df.currentRowIdx);
                let currentPatientVisit = this.totalVisits[df.currentCol.name][subjId]
                currentPatientVisit.updateSubjectVisitDomains();
                grok.shell.o = await this.patientVisitPanel(currentPatientVisit);
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
            if (study.domains[domain] && [SUBJECT_ID, VISIT_NAME, this.proceduresAtVisit[domain].column].every(it => study.domains[domain].columns.names().includes(it))) {
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
                    this.totalVisits[colName][subjId] = new PatientVisit(study.domains);
                    this.totalVisits[colName][subjId].updateSubjectVisit(subjId, currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
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
                    if (!study.domains[domain].col(this.eventsSinceLastVisit[domain].column).isNone(i)) {
                        let startDay = study.domains[domain].get(this.eventsSinceLastVisit[domain].column, i);
                        let subjId = study.domains[domain].get(SUBJECT_ID, i);
                        this.subjSet.add(subjId);
                        for (let z = 0; z < this.sortedVisitNamesAndDays.length - 1; z++) {
                            if (startDay > this.sortedVisitNamesAndDays[z].day && startDay <= this.sortedVisitNamesAndDays[z + 1].day) {
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

    async studyVisitPanel() {
        let panelDiv = ui.div();
        let acc = this.createAccWithTitle('Study visit panel', `${this.studyVisit.currentVisitName}`);

        acc.addPane(`Summary`, () => {
            return ui.tableFromMap({
                'Visit day': this.studyVisit.currentVisitDay,
                'Total patients': this.studyVisit.totalPatients,
                'Min visit date': this.studyVisit.minVisitDate,
                'Max visit dae': this.studyVisit.maxVisitDate,
            })
        });

        let getRowNumber = (df) => {
            return df ? df.rowCount === 1 && df.getCol(SUBJECT_ID).isNone(0) ? 0 : df.rowCount : 0;
        }

        if (this.studyVisit.exAtVisit && this.studyVisit.exAtVisit.columns.names().includes(this.studyVisit.exAtVisit)) {
            acc.addPane('Drug exposure', () => {
                if (!getRowNumber(this.studyVisit.exAtVisit)) {
                    return ui.divText('No records found');
                }
                return DG.Viewer.fromType(DG.VIEWER.PIE_CHART, this.studyVisit.exAtVisit, {
                    category: this.studyVisit.extrtWithDoseColName,
                }).root;
            })
        }

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

        let aeRowNum = getRowNumber(this.studyVisit.aeSincePreviusVisit);
        createPane('AEs since last visit', aeRowNum, this.studyVisit.aeSincePreviusVisit, AE_DECOD_TERM);

        let cmRowNum = getRowNumber(this.studyVisit.conmedSincePreviusVisit);
        createPane('CONMEDs since last visit', cmRowNum, this.studyVisit.conmedSincePreviusVisit, CON_MED_NAME);

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
                            showValueSelector: false
                        });
                        return plot.root;
                    })
                })
                return categoriesAcc.root;
            })
        }
        if (this.studyVisit.lbAtVisit && [LAB_TEST, LAB_RES_N].every(it => this.studyVisit.lbAtVisit.columns.names().includes(it))) {
            createDistributionPane('Laboratory', this.studyVisit.lbAtVisit, LAB_TEST, LAB_RES_N);
        }
        if (this.studyVisit.vsAtVisit && [VS_TEST, VS_RES_N].every(it => this.studyVisit.vsAtVisit.columns.names().includes(it))) {
            createDistributionPane('Vital signs', this.studyVisit.vsAtVisit, VS_TEST, VS_RES_N);
        }

        panelDiv.append(acc.root);

        return panelDiv;

    }

    async patientVisitPanel(patientVisit: PatientVisit) {
        let panelDiv = ui.div();
        let acc = this.createAccWithTitle('Patient visit panel', `${patientVisit.currentSubjId}`);
        
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
                if (patientVisit[domain]) {
                    const rowNum =
                        patientVisit[domain].rowCount === 1 && patientVisit[domain].getCol(SUBJECT_ID).isNone(0) ?
                            0 : patientVisit[domain].rowCount
                    createPane(domain, key, rowNum);
                }
            });
        }

        createAccordion();

        panelDiv.append(acc.root);

        return panelDiv;

    }

}
