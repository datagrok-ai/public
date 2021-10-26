import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain, getMaxVisitName, getMinVisitName, getUniqueValues, getVisitNamesAndDays } from '../data-preparation/utils';
import { ETHNIC, RACE, SEX, TREATMENT_ARM } from '../constants';
import { labDynamicComparedToBaseline } from '../data-preparation/data-preparation';
import { ILazyLoading } from '../lazy-loading/lazy-loading';


export class TimeProfileView extends DG.ViewBase implements ILazyLoading {

    blVisitChoices: DG.InputBase;
    epVisitChoices: DG.InputBase;
    laboratoryDataFrame: DG.DataFrame;
    relativeChangeFromBlDataFrame: DG.DataFrame;
    uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, 'LBTEST'));
    uniqueVisits = Array.from(getUniqueValues(study.domains.lb, 'VISIT'));
    splitBy = [ TREATMENT_ARM, SEX, RACE, ETHNIC ];
    types = ['Values', 'Changes'];
    selectedLabValue: string;
    selectedType: string;
    bl: string;
    ep: string;
    visitNamesAndDays: any [];
    linechart: any;

    constructor(name) {
        super({});
        this.name = name;
        this.helpUrl = 'https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/views_help/time_profile.md';
    }

    loaded: boolean;

    load(): void {
        this.selectedLabValue = this.uniqueLabValues[ 0 ] as string;
        this.selectedType = this.types[0];
        this.visitNamesAndDays = getVisitNamesAndDays(study.domains.lb);
        this.bl = this.visitNamesAndDays[0].name;
        this.ep = this.visitNamesAndDays[this.visitNamesAndDays.length-1].name;
        this.createLaboratoryDataframe();
        this.createrelativeChangeFromBlDataframe();

/*         let labValueChoices = ui.choiceInput('BL', this.selectedLabValue, this.uniqueLabValues);
        labValueChoices.onChanged((v) => {
            this.selectedLabValue = labValueChoices.value;
            this.updateTimeProfileDiv();
        }); */

        let typeChoices = ui.choiceInput('', this.selectedType, this.types);
        typeChoices.onChanged((v) => {
            this.selectedType = typeChoices.value;
            this.updateTimeProfile();
        });

        this.blVisitChoices = ui.choiceInput('', this.bl, this.uniqueVisits);
        this.blVisitChoices.onChanged((v) => {
            this.bl = this.blVisitChoices.value;
            this.updateTimeProfile();
        });

        this.epVisitChoices = ui.choiceInput('', this.ep, this.uniqueVisits);
        this.epVisitChoices.onChanged((v) => {
            this.ep = this.epVisitChoices.value;
            this.updateTimeProfile();
        });

        this.root.className = 'grok-view ui-box';
/*         this.root.append(ui.splitV([
            ui.box(ui.divH([ //labValueChoices.root, 
                typeChoices.root, 
                this.blVisitChicesDiv
            ]), { style: { maxHeight: '100px' } }),
            this.timeProfileDiv
        ])) */
        this.linechart = DG.Viewer.lineChart(this.filterDataFrameByDays(this.laboratoryDataFrame), {
            splitColumnName: this.splitBy[0],
            xColumnName: 'VISITDY',
            yColumnNames: [`${this.selectedLabValue} avg(LBSTRESN)`],
            whiskersType: 'Med | Q1, Q3'
        });
        this.root.append(this.linechart.root);
        this.setRibbonPanels([
            [
                ui.span([ 'Plot ' ]),
                typeChoices.root,
                ui.span([' from ']),
                this.blVisitChoices.root,
                ui.span([' to ']),
                this.epVisitChoices.root
            ]
        ]);
    }

    private updateTimeProfile() {
        switch (this.selectedType) {
            case 'Values': {
                this.linechart.dataFrame = this.filterDataFrameByDays(this.laboratoryDataFrame);
                break;
            }
            case 'Changes': {
                this.createrelativeChangeFromBlDataframe();
                this.linechart.dataFrame = this.relativeChangeFromBlDataFrame;
                break;
            }
            default: {
                break;
            }
        }
    }

    private createLaboratoryDataframe() {
        let dfWithArm = addDataFromDmDomain(study.domains.lb, study.domains.dm, [ 'USUBJID', 'VISITDY', 'VISIT', 'LBTEST', 'LBSTRESN' ], this.splitBy);
        this.laboratoryDataFrame = this.createPivotedDataframe(dfWithArm, 'LBSTRESN');
    }

    private createrelativeChangeFromBlDataframe(){
        let df = this.filterDataFrameByDays(study.domains.lb.clone());
        labDynamicComparedToBaseline(df,  this.bl, 'VISIT', 'LAB_DYNAMIC_BL', true);
        let dfWithArm = addDataFromDmDomain(df, study.domains.dm, [ 'USUBJID', 'VISITDY', 'VISIT', 'LBTEST', 'LBSTRESN' ], this.splitBy);
        this.relativeChangeFromBlDataFrame = this.createPivotedDataframe(dfWithArm, 'LBSTRESN');
    }

    private createPivotedDataframe(df: DG.DataFrame, aggregatedColName: string) {
        return df
            .groupBy([ 'USUBJID', 'VISITDY' ].concat(this.splitBy))
            .pivot('LBTEST')
            .avg(aggregatedColName)
            .aggregate();
    }

    private filterDataFrameByDays(df: DG.DataFrame){
        let blDay = this.visitNamesAndDays.find(it => it.name === this.bl).day;
        let epDay = this.visitNamesAndDays.find(it => it.name === this.ep).day;
        let filteredDf = df.groupBy(df.columns.names())
        .where(`VISITDY >= ${blDay} and VISITDY <= ${epDay}`)
        .aggregate();
        return filteredDf;
    }

}