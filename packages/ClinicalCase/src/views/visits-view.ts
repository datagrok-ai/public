import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { updateDivInnerHTML } from '../utils/utils';
import { createPivotedDataframe, getVisitNamesAndDays, addDataFromDmDomain } from '../data-preparation/utils';
import { LAB_RES_N, LAB_TEST, SUBJECT_ID, VISIT_DAY, VISIT_NAME, VISIT_START_DATE, VS_RES_N, VS_TEST } from '../constants/columns-constants';
import { PatientVisit } from '../model/patient-visit';
import { StudyVisit } from '../model/study-visit';
import { _package } from '../package';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';
import { AE_START_DAY_FIELD, AE_TERM_FIELD, CON_MED_NAME_FIELD, CON_MED_START_DAY_FIELD, INV_DRUG_NAME_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from '../views-config';
import { VISITS_VIEW_NAME } from '../constants/view-names-constants';
import { DOMAINS_COLOR_PALETTE } from '../constants/constants';

export class VisitsView extends ClinicalCaseViewBase {

    sv: DG.DataFrame;
    pivotedSv: DG.DataFrame;
    sortedVisitNamesAndDays: any;
    sortedVisitNames: any;
    patientVisit = new PatientVisit();
    studyVisit = new StudyVisit();
    totalVisits = {};
    proceduresAtVisit: any;
    eventsSinceLastVisit: any;
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
    switchGrid: any;

    constructor(name) {
        super({});
        this.name = name;
        this.helpUrl = `${_package.webRoot}/views_help/visits.md`;
    }

    createView(): void {
        this.proceduresAtVisit = this.getProceduresAtVisitDict();
        this.eventsSinceLastVisit = this.eventsSinceLastVisitCols();
        this.existingDomains = Object.keys(this.proceduresAtVisit)
            .concat(Object.keys(this.eventsSinceLastVisit))
            .filter(it => study.domains[it] !== null);
        this.assignColorsToDomains();
        this.sortedVisitNamesAndDays = getVisitNamesAndDays(study.domains.sv, VISIT_NAME, VISIT_DAY, true);
        this.sortedVisitNames = this.sortedVisitNamesAndDays.map(it => it.name);

        this.switchGrid = this.createSwitchGridInput();
        this.createHeatMapRibbon();
        this.createGridRibbon();
        updateDivInnerHTML(this.ribbonDiv, this.gridRibbonDiv);

        this.updateGridAndHeatMap();

        this.root.className = 'grok-view ui-box';
        this.root.append(this.div);

        this.setRibbonPanels(
            [[this.switchGrid.root], [this.ribbonDiv], [this.visitsConfig()]] ,
        );
    }

    private updateGridAndHeatMap() {
        this.createPivotedSv();
        this.createTotalVisits();
        this.visitsGrid = this.pivotedSv.plot.grid();
        this.renderVisitCell();
        this.heatMap = DG.Viewer.fromType(DG.VIEWER.HEAT_MAP, this.createHeatMapDf(), {
            'colorCoding': 'All',
        });
        this.switchGrid.value = true;
        updateDivInnerHTML(this.div, this.visitsGrid.root);
    }

    private createPivotedSv() {
        this.sv = study.domains.sv.clone();
        this.pivotedSv = createPivotedDataframe(this.sv, [SUBJECT_ID], VISIT_NAME, VISIT_START_DATE, []);
        this.pivotedSv.columns.names().forEach(col => {
            if (this.pivotedSv.getCol(col).name !== VISIT_START_DATE) {
                this.pivotedSv.getCol(col).tags.format = 'yyyy-MM-dd';
            }
            this.pivotedSv.getCol(col).name = col.replace(` first(${VISIT_START_DATE})`, '');
        });
        this.pivotedSv = addDataFromDmDomain(this.pivotedSv, study.domains.dm, this.pivotedSv.columns.names(), [VIEWS_CONFIG[this.name][TRT_ARM_FIELD]]);
        this.pivotedSv = this.pivotedSv.clone(null, [SUBJECT_ID, VIEWS_CONFIG[this.name][TRT_ARM_FIELD]].concat(this.sortedVisitNames));
        this.pivotedSv.onCurrentCellChanged.subscribe(() => {
            setTimeout(() => {
                this.createVisitPropertyPanel(this.pivotedSv);

            }, 100);

        });
    }

    private eventsSinceLastVisitCols() {
        const eventsSinceLastVisit = { 'ae': { column: VIEWS_CONFIG[this.name][AE_START_DAY_FIELD] }, 'cm': { column: VIEWS_CONFIG[this.name][CON_MED_START_DAY_FIELD] } };
        const filtered = Object.assign({}, ...
            Object.entries(eventsSinceLastVisit).filter(([k,v]) => study.domains[k].col(v.column)).map(([k,v]) => ({[k]:v}))
        );
        return filtered;
    }

    private visitsConfig() {
        let visitsConfig = ui.iconFA('cog', () => {
            let div = ui.div();
            let visitsDf = DG.DataFrame.create();
            visitsDf.columns.addNewString('Name');
            visitsDf.columns.addNewInt('Day');
            visitsDf.columns.addNewString('Action');
            this.sortedVisitNamesAndDays.forEach(visit => {
                visitsDf.rows.addNew([visit.name, visit.day, ''], false);
            });
            let grid = visitsDf.plot.grid();
            grid.setOptions({ allowRowSelection: false, extendLastColumn: true });
            let col = grid.columns.byName('Action');
            col.cellType = 'html';

            grid.onCellPrepare(function (gc) {
                if (gc.isTableCell && gc.tableColumn.name === 'Action') {
                    let eventElement = ui.icons.delete(() => {}, 'Delete');
                    gc.style.element = ui.button(eventElement, () => {
                        gc.grid.dataFrame.rows.removeAt(gc.gridRow);
                        let scrollTo = gc.grid.dataFrame.rowCount === gc.gridRow ? gc.gridRow - 1 : gc.gridRow;
                        grid.scrollToCell('Action', scrollTo);
                    });
                    gc.style.element.style.paddingBottom = '7px';
                    gc.style.element.style.paddingLeft = '15px';
                }
            });

            let addButton = ui.button(ui.icons.add(() => { }), () => {
                visitsDf.rows.addNew();
            });
            div.append(ui.div(grid.root));
            div.append(addButton);
            ui.dialog({ title: 'Visits' })
                .add(div)
                .onOK(() => {
                    this.totalVisits = {};
                    this.sortedVisitNamesAndDays = getVisitNamesAndDays(visitsDf, 'Name', 'Day', true);
                    this.sortedVisitNames = this.sortedVisitNamesAndDays.map(it => it.name);
                    this.updateGridAndHeatMap();
                })
                .show();
        });
        return visitsConfig;
    }

    private getProceduresAtVisitDict() {
        return { 'lb': { column: LAB_TEST }, 'ex': { column: VIEWS_CONFIG[VISITS_VIEW_NAME][INV_DRUG_NAME_FIELD] }, 'vs': { column: VS_TEST } };
    }

    private assignColorsToDomains() {
        this.existingDomains.forEach((domain) => {
            this.proceduresAtVisit[domain] ?
                this.proceduresAtVisit[domain]['color'] = DG.Color.toRgb(DOMAINS_COLOR_PALETTE[domain]) :
                this.eventsSinceLastVisit[domain]['color'] = DG.Color.toRgb(DOMAINS_COLOR_PALETTE[domain]);
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
        this.selectedDomainsDiv.addEventListener('click', (event) => {
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
            event.stopPropagation();
          });

        this.gridRibbonDiv = ui.divH([
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
         //   domainName.style.paddingTop = '10px';
            selectedDomainsDiv.append(domainName);
        });
        return selectedDomainsDiv;
    }


    private createHeatMapDf() {
        let df = DG.DataFrame.create();
        typeof (study.domains.dm.get(SUBJECT_ID, 0)) === 'number' ? df.columns.addNewInt(SUBJECT_ID) : df.columns.addNewString(SUBJECT_ID);
        this.pivotedSv.columns.names().forEach(col => {
            if (col !== SUBJECT_ID && col !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]) {
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
        df = addDataFromDmDomain(df, study.domains.dm, df.columns.names().filter(it => it !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]), [VIEWS_CONFIG[this.name][TRT_ARM_FIELD]]);
        df = df.clone(null, [SUBJECT_ID, VIEWS_CONFIG[this.name][TRT_ARM_FIELD]].concat(this.sortedVisitNames));
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
            if (col !== SUBJECT_ID && col !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]) {
                df.col(col).tags[DG.TAGS.COLOR_CODING_TYPE] = 'Linear';
                df.col(col).tags[DG.TAGS.COLOR_CODING_LINEAR] = `[${DG.Color.white}, ${DG.Color.blue}]`;
            }
        });
    }

    private async createVisitPropertyPanel(df: DG.DataFrame) {
        if (df.currentCol.name !== SUBJECT_ID && df.currentCol.name !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]) {
            if (df.currentRowIdx === -1) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(df.currentCol.name);
                this.studyVisit.updateStudyVisit(study.domains, currentVisit.day, currentVisit.name, previousVisit ? previousVisit.day : null);
                grok.shell.o = await this.studyVisitPanel();
            } else {
                let subjId = df.get(SUBJECT_ID, df.currentRowIdx);
                let currentPatientVisit = this.totalVisits[df.currentCol.name][subjId]
                currentPatientVisit.updateSubjectVisitDomains(study.domains);
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
            if (colName !== SUBJECT_ID && colName !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]) {
                let { current: currentVisit, previous: previousVisit } = this.getCurrentAndPreviousVisits(colName);
                this.totalVisits[colName] = {};

                this.pivotedSv.getCol(SUBJECT_ID).categories.forEach(subjId => {
                    this.totalVisits[colName][subjId] = new PatientVisit();
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
            if (gc.isTableCell && gc.gridColumn.name !== SUBJECT_ID && gc.gridColumn.name !== VIEWS_CONFIG[this.name][TRT_ARM_FIELD]) {
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
        let acc = this.createAccWithTitle('Study visit panel', `${this.studyVisit.name}`);

        acc.addPane(`Summary`, () => {
            return ui.tableFromMap({
                'Visit day': this.studyVisit.day,
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

        let createSinceLastVisitPane = (df, col, paneName) => {
            if(df) {
                let aeRowNum = getRowNumber(df);
                createPane(paneName, aeRowNum, df, col);
            }
        }

        createSinceLastVisitPane(this.studyVisit.aeSincePreviusVisit, VIEWS_CONFIG[this.name][AE_TERM_FIELD], 'AEs since last visit');
        createSinceLastVisitPane(this.studyVisit.conmedSincePreviusVisit, VIEWS_CONFIG[this.name][CON_MED_NAME_FIELD], 'CONMEDs since last visit');

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
