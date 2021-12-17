import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { requiredColumnsByView } from '../constants';
import { MH_BODY_SYSTEM, MH_CATEGORY, MH_DECOD_TERM, MH_TERM, SUBJECT_ID } from '../columns-constants';

export class MedicalHistoryView extends DG.ViewBase implements ILazyLoading {

    mh: DG.DataFrame;
    mhReplacedTermColName = 'mh_term';

    constructor(name) {
        super({});
        this.name = name;
    }

    loaded: boolean;

    load(): void {
        checkMissingDomains(requiredColumnsByView[this.name], this);
    }

    createView(): void {

        this.mh = study.domains.mh.clone();
        this.replaceNaDecodeTerm();

        let viewerTitle = {
            style: {
                'color': 'var(--grey-6)',
                'margin': '12px 0px 6px 12px',
                'font-size': '16px',
            }
        };


        let grid = this.mh.plot.grid();

        let mhCategoryPie = DG.Viewer.fromType(DG.VIEWER.PIE_CHART, this.mh, {
            category: MH_CATEGORY,
        });
        mhCategoryPie.root.prepend(ui.divText('By category', viewerTitle));

        let mhDecodTermChart = this.mh.plot.bar({
            split: this.mhReplacedTermColName,
            style: 'dashboard',
            legendVisibility: 'Never'
        }).root;
        mhDecodTermChart.prepend(ui.divText('By term', viewerTitle));

        let mhBodySystemChart = this.mh.plot.bar({
            split: MH_BODY_SYSTEM,
            style: 'dashboard',
            legendVisibility: 'Never'
        }).root;
        mhBodySystemChart.prepend(ui.divText('By body system', viewerTitle));

        this.root.className = 'grok-view ui-box';
        this.root.append(ui.splitV([
            ui.splitH([
                mhDecodTermChart,
                mhBodySystemChart
            ]),
            ui.splitH([
                mhCategoryPie.root,
                grid.root
            ])
        ]))

    }

    private replaceNaDecodeTerm() {
        this.mh.columns.addNewString(this.mhReplacedTermColName)
            .init((i) => this.mh.getCol(MH_DECOD_TERM).isNone(i) ? this.mh.get(MH_TERM, i) : this.mh.get(MH_DECOD_TERM, i));
    }
}