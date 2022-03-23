import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import { study } from '../clinical-study';
import { AE_END_DATE, AE_SEVERITY, AE_START_DATE, AGE, RACE, SEX, SUBJECT_ID } from "../constants/columns-constants";
import * as ui from 'datagrok-api/ui';
import { dictToString, getNullOrValue } from "../data-preparation/utils";
import { getSubjectDmData } from "../data-preparation/data-preparation";
import { SEVERITY_COLOR_DICT } from "../constants/constants";
import { updateDivInnerHTML } from "../utils/utils";
import { AE_END_DAY_FIELD, AE_START_DAY_FIELD, AE_TERM_FIELD, TRT_ARM_FIELD, VIEWS_CONFIG } from "../views-config";
import { AE_BROWSER_VIEW_NAME } from "../constants/view-names-constants";

export class AEBrowserHelper {

    domains = ['ae', 'ex', 'cm'];
    dyDomains = ['vs', 'lb'];
    domainsWithoutVisitDays = ['mh'];
    domainsToExclude = ['dm', 'tv', 'sv'];
    additionalDomains = [];
    selectedAdditionalDomains = [];
    aeToSelect: DG.DataFrame;
    daysPriorAe = 5;
    currentSubjId = '';
    currentAeDay: number;
    name = AE_BROWSER_VIEW_NAME;

    constructor(dataFrame: DG.DataFrame) {
        let presentDomains = [];
        this.aeToSelect = dataFrame;
        study.domains.all().forEach(it => {
            if (!this.domainsToExclude.includes(it.name)) {
                this[it.name] = study.domains[it.name].clone();
                presentDomains.push(it.name);
                if (!this.domains.includes(it.name)) {
                    this.additionalDomains.push(it.name);
                }
            }
        });
        this.domains = this.domains.filter(it => presentDomains.includes(it));
    }

    updateDomains() {
        this.domains.concat(this.selectedAdditionalDomains).forEach(domain => {
            const condition = this.domainsWithoutVisitDays.includes(domain) ?
                `${SUBJECT_ID} = ${this.currentSubjId}` :
                this.dyDomains.includes(domain) ?
                    study.domains[domain].col(`${domain.toUpperCase()}DY`) ?
                        `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}DY < ${this.currentAeDay} and ${domain.toUpperCase()}DY > ${this.currentAeDay - this.daysPriorAe}` :
                        null :
                    [`${domain.toUpperCase()}STDY`, `${domain.toUpperCase()}ENDY`].every(it => study.domains[domain].col(it) !== null) ?
                        `${SUBJECT_ID} = ${this.currentSubjId} and ${domain.toUpperCase()}STDY < ${this.currentAeDay} and ${domain.toUpperCase()}ENDY > ${this.currentAeDay - this.daysPriorAe}` :
                        null;
            if (condition) {
                this[domain] = study.domains[domain]
                .clone()
                .groupBy(study.domains[domain].columns.names())
                .where(`${condition}`)
                .aggregate();
            } else {
                this[domain] = null;
            }
        })
    }

    async propertyPanel() {
        this.updateDomains();
        grok.shell.o = await this.aeBrowserPanel();
    }

    async aeBrowserPanel() {

        let panelDiv = ui.div();
        let accae = ui.accordion(`${this.name} panel`);
        let accIcon = ui.element('i');
        accIcon.className = 'grok-icon svg-icon svg-view-layout';

        if (this.aeToSelect.currentRowIdx !== -1) {
            let subjId = this.aeToSelect.get(SUBJECT_ID, this.aeToSelect.currentRowIdx);
            let title = ui.tooltip.bind(ui.label(subjId), dictToString(getSubjectDmData(subjId, [AGE, SEX, RACE, VIEWS_CONFIG[this.name][TRT_ARM_FIELD]])));
            let description = ui.divH([ui.divText(String(this.aeToSelect.get(VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_TERM_FIELD], this.aeToSelect.currentRowIdx).toLowerCase()))]);

            if (this.aeToSelect.columns.names().includes(AE_SEVERITY)) {
                let severity = this.aeToSelect.get(AE_SEVERITY, this.aeToSelect.currentRowIdx);
                let severityStyle = {
                    style: {
                        color: `${SEVERITY_COLOR_DICT[severity.toUpperCase()]}`,
                        marginRight: '5px',
                        fontWeight: 'bold'
                    }
                }
                description.prepend(ui.divText(severity, severityStyle));
            }

            accae.addTitle(ui.span([accIcon, title]));

            let getAeDate = (col) => {
                this.aeToSelect.col(col) ? getNullOrValue(this.aeToSelect, AE_START_DATE, this.aeToSelect.currentRowIdx) : 'null';
            }

            let startEndDays = ui.tooltip.bind(
                ui.label(`${getNullOrValue(this.aeToSelect, VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_START_DAY_FIELD], this.aeToSelect.currentRowIdx)} - ${getNullOrValue(this.aeToSelect, VIEWS_CONFIG[AE_BROWSER_VIEW_NAME][AE_END_DAY_FIELD], this.aeToSelect.currentRowIdx)}`),
                `${getAeDate(AE_START_DATE)} - ${getAeDate(AE_END_DATE)}`);


            let daysInput = ui.intInput('Prior AE', this.daysPriorAe);
            daysInput.onChanged((v) => {
                this.daysPriorAe = daysInput.value;
                this.updateDomains();
                updateAccordion();
            });
            startEndDays.innerHTML = 'Days ' + startEndDays.innerHTML;
            startEndDays.style.marginTop = '5px';
            //@ts-ignore
            accae.header = ui.div([
                description,
                startEndDays,
                //@ts-ignore
                ui.divH([ui.divText('Days prior AE'), daysInput.input], { style: { alignItems: 'center', gap: '5px' } }),
            ])

            let getPaneContent = (it, rowNum) => {
                if (it) {
                    if (!rowNum) {
                        return ui.divText('No records found');
                    } else {
                        let grid = this[it].plot.grid();
                        if (rowNum < 7) {
                            grid.root.style.maxHeight = rowNum < 4 ? '100px' : '150px';
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
                const totalDomains = this.domains.concat(this.selectedAdditionalDomains);
                const panesToRemove = accae.panes.filter(it => !totalDomains.includes(it.name));
                panesToRemove.forEach(it => accae.removePane(it));
                totalDomains.forEach(it => {
                    if (this[it]) {
                        const rowNum = this[it].rowCount === 1 && this[it].getCol(SUBJECT_ID).isNone(0) ? 0 : this[it].rowCount
                        let pane = accae.getPane(`${it}`);
                        if (pane) {
                            //@ts-ignore
                            updateDivInnerHTML(pane.root.lastChild, getPaneContent(it, rowNum));
                            //@ts-ignore
                            pane.root.firstChild.lastChild.innerText = rowNum;
                        } else {
                            createPane(it, rowNum);
                        }
                    }
                })
            }

            let createAccordion = () => {
                this.domains.concat(this.selectedAdditionalDomains).forEach(it => {
                    if (this[it]) {
                        const rowNum = this[it].rowCount === 1 && this[it].getCol(SUBJECT_ID).isNone(0) ? 0 : this[it].rowCount;
                        createPane(it, rowNum);
                    }
                })
            }

            let addButton = ui.button(ui.icons.add(() => { }), () => {
                let domainsMultiChoices = ui.multiChoiceInput('', this.selectedAdditionalDomains, this.additionalDomains);
                ui.dialog({ title: 'Select domains' })
                    .add(ui.div([domainsMultiChoices]))
                    .onOK(() => {
                        this.selectedAdditionalDomains = domainsMultiChoices.value;
                        this.updateDomains();
                        updateAccordion();
                    })
                    //@ts-ignore
                    .show({ centerAt: addButton });
            });
            createAccordion();

            panelDiv.append(accae.root);
            panelDiv.append(addButton);

            return panelDiv;
        }
        return accae.root;
    }

}