import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { UIBuilder } from '../helpers/ui-builder';

export class ConditionOccurrenceAchillesView extends DG.ViewBase {

    constructor(name) {
        super({});
        this.name = name;
        this.createView();
    }

    async createView(): Promise<void> {
        let uiBuilder = new UIBuilder(this.name, this.root, [ 'PrevalenceByGenderAgeYear', 'PrevalenceByMonth', 'ByType' ]);
    }
}