import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study, ClinRow } from "../clinical-study";
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { requiredColumnsByView } from '../constants';
import { addDataFromDmDomain } from '../data-preparation/utils';
import { ETHNIC, RACE, SEX, TREATMENT_ARM } from '../columns-constants';

export class TreeMapView extends DG.ViewBase implements ILazyLoading {

    aeDataframeWithDm: DG.DataFrame
    dmFields = [TREATMENT_ARM, SEX, RACE, ETHNIC];
    selectedSplitBy = '';
    treeMap: any;

    constructor(name) {
        super({});
        this.name = name;
    }

    loaded: boolean;

    load(): void {
        checkMissingDomains(requiredColumnsByView[this.name], this);
    }

    createView(): void {
        this.dmFields = this.dmFields.filter(it => study.domains.dm.columns.names().includes(it));
        this.aeDataframeWithDm = addDataFromDmDomain(study.domains.ae, study.domains.dm, study.domains.ae.columns.names(), this.dmFields);

        this.treeMap = DG.Viewer.fromType(DG.VIEWER.TREE_MAP, this.aeDataframeWithDm, {
            "splitByColumnNames": [
                this.dmFields[0],
                ""
            ],
            "colorAggrType": "count"
          });
          this.root.className = 'grok-view ui-box';
          this.root.appendChild(this.treeMap.root);

    }
}