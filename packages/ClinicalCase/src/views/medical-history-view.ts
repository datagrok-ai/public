import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { MH_BODY_SYSTEM, MH_CATEGORY, MH_DECOD_TERM, MH_TERM } from '../columns-constants';
import { _package } from '../package';
import { ClinicalCaseViewBase } from '../model/ClinicalCaseViewBase';
import { checkColumnsAndCreateViewer, updateDivInnerHTML } from "./utils";

export class MedicalHistoryView extends ClinicalCaseViewBase {

    mh: DG.DataFrame;
    mhReplacedTermColName = 'mh_term';
    mhCategoryPie = ui.box();
    mhDecodTermChart = ui.box();
    mhBodySystemChart = ui.box();

    constructor(name) {
        super({});
        this.name = name;
        this.helpUrl = `${_package.webRoot}/views_help/medical_history.md`;
    }

    createView(): void {

        this.mh = study.domains.mh.clone();
        if (this.mh.col(MH_DECOD_TERM) && this.mh.col(MH_TERM)) {
            this.replaceNaDecodeTermWihTerm();
        } else {
            this.mhReplacedTermColName = MH_DECOD_TERM;
        }

        let viewerTitle = {
            style: {
                'color': 'var(--grey-6)',
                'margin': '12px 0px 6px 12px',
                'font-size': '16px',
            }
        };

        let grid = this.mh.plot.grid();

        checkColumnsAndCreateViewer(
            study.domains.mh,
            [MH_CATEGORY],
            this.mhCategoryPie, () => {
                this.createMhCategoryPie(viewerTitle);
            },
            'By category');

        checkColumnsAndCreateViewer(
            study.domains.mh,
            [MH_DECOD_TERM],
            this.mhDecodTermChart, () => {
                this.createDecodeTermChart(viewerTitle);
            },
            'By term');

        checkColumnsAndCreateViewer(
            study.domains.mh,
            [MH_BODY_SYSTEM],
            this.mhBodySystemChart, () => {
                this.createBodySystemChart(viewerTitle);
            },
            'By body system');

        this.root.className = 'grok-view ui-box';
        this.root.append(ui.splitV([
            ui.splitH([
                this.mhDecodTermChart,
                this.mhBodySystemChart
            ]),
            ui.splitH([
                this.mhCategoryPie,
                grid.root
            ])
        ]))

    }

    private replaceNaDecodeTermWihTerm() {
        this.mh.columns.addNewString(this.mhReplacedTermColName)
            .init((i) => this.mh.getCol(MH_DECOD_TERM).isNone(i) ? this.mh.get(MH_TERM, i) : this.mh.get(MH_DECOD_TERM, i));
    }

    private createMhCategoryPie(viewerTitle: any) {

        let mhCategoryPie = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, this.mh, {
            category: MH_CATEGORY,
        });
        mhCategoryPie.root.prepend(ui.divText('By category', viewerTitle));
        updateDivInnerHTML(this.mhCategoryPie, mhCategoryPie.root);
    }

    private createDecodeTermChart(viewerTitle: any) {
        let mhDecodTermChart = this.mh.plot.bar({
            split: this.mhReplacedTermColName,
            style: 'dashboard',
            legendVisibility: 'Never'
        }).root;
        mhDecodTermChart.prepend(ui.divText('By term', viewerTitle));
        updateDivInnerHTML(this.mhDecodTermChart, mhDecodTermChart);
    }

    private createBodySystemChart(viewerTitle: any) {
        let mhBodySystemChart = this.mh.plot.bar({
            split: MH_BODY_SYSTEM,
            style: 'dashboard',
            legendVisibility: 'Never'
        }).root;
        mhBodySystemChart.prepend(ui.divText('By body system', viewerTitle));
        updateDivInnerHTML(this.mhBodySystemChart, mhBodySystemChart);
    }
}