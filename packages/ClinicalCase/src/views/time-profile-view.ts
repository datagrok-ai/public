import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addDataFromDmDomain, getUniqueValues, getLabVisitNamesAndDays } from '../data-preparation/utils';
import { ETHNIC, LAB_RES_N, LAB_TEST, LAB_VISIT_DAY, LAB_VISIT_NAME, RACE, SEX, SUBJECT_ID, TREATMENT_ARM } from '../columns-constants';
import { labDynamicComparedToBaseline } from '../data-preparation/data-preparation';
import { ILazyLoading } from '../lazy-loading/lazy-loading';
import { checkMissingDomains } from './utils';
import { _package } from '../package';
import { requiredColumnsByView } from '../constants';


export class TimeProfileView extends DG.ViewBase implements ILazyLoading {

    blVisitChoices: DG.InputBase;
    epVisitChoices: DG.InputBase;
    laboratoryDataFrame: DG.DataFrame;
    relativeChangeFromBlDataFrame: DG.DataFrame;
    uniqueLabValues: any;
    uniqueVisits: any;
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
        this.helpUrl = `${_package.webRoot}/views_help/time_profile.md`;
    }

    loaded: boolean;

    load(): void {
        checkMissingDomains(requiredColumnsByView[this.name], false, this);
     }

    createView(): void {
        this.uniqueLabValues = Array.from(getUniqueValues(study.domains.lb, LAB_TEST));
        this.uniqueVisits = Array.from(getUniqueValues(study.domains.lb, LAB_VISIT_NAME));
        this.selectedLabValue = this.uniqueLabValues[ 0 ] as string;
        this.selectedType = this.types[0];
        this.visitNamesAndDays = getLabVisitNamesAndDays(study.domains.lb);
        this.bl = this.visitNamesAndDays[0].name;
        this.ep = this.visitNamesAndDays[this.visitNamesAndDays.length-1].name;
        this.createLaboratoryDataframe();

        let typeChoices = ui.choiceInput('', this.selectedType, this.types);
        typeChoices.onChanged((v) => {
            this.selectedType = typeChoices.value;
            this.updateTimeProfile();
        });

        let labChoices = ui.choiceInput('', this.selectedLabValue, this.uniqueLabValues);
        labChoices.onChanged((v) => {
            this.selectedLabValue = labChoices.value;
            this.updateTimeProfile();
        });
        //@ts-ignore
        labChoices.input.style.width = '200px';

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
        this.linechart = DG.Viewer.lineChart(this.laboratoryDataFrame, {
            splitColumnName: this.splitBy[0],
            xColumnName: LAB_VISIT_DAY,
            yColumnNames: [`${this.selectedLabValue} avg(${LAB_RES_N})`],
            whiskersType: 'Med | Q1, Q3'
        });
        this.root.append(this.linechart.root);
        this.setRibbonPanels([
            [
                ui.span([ 'Plot ' ]),
                labChoices.root,
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
                this.createLaboratoryDataframe();
                this.linechart.dataFrame = this.laboratoryDataFrame;
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
        let df = this.filterDataFrameByDays(study.domains.lb.clone());
        let dfWithArm = addDataFromDmDomain(df, study.domains.dm, [ SUBJECT_ID, LAB_VISIT_DAY, LAB_VISIT_NAME, LAB_TEST, LAB_RES_N ], this.splitBy);
        this.laboratoryDataFrame = this.createPivotedDataframe(dfWithArm, LAB_RES_N);
    }

    private createrelativeChangeFromBlDataframe(){
        let df = this.filterDataFrameByDays(study.domains.lb.clone());
        labDynamicComparedToBaseline(df,  this.bl, LAB_VISIT_NAME, 'LAB_DYNAMIC_BL', true);
        let dfWithArm = addDataFromDmDomain(df, study.domains.dm, [ SUBJECT_ID, LAB_VISIT_DAY, LAB_VISIT_NAME, LAB_TEST, LAB_RES_N ], this.splitBy);
        this.relativeChangeFromBlDataFrame = this.createPivotedDataframe(dfWithArm, LAB_RES_N);
    }

    private createPivotedDataframe(df: DG.DataFrame, aggregatedColName: string) {
        return df
            .groupBy([ SUBJECT_ID, LAB_VISIT_DAY ].concat(this.splitBy))
            .pivot(LAB_TEST)
            .avg(aggregatedColName)
            .aggregate();
    }

    private filterDataFrameByDays(df: DG.DataFrame){
        let blDay = this.visitNamesAndDays.find(it => it.name === this.bl).day;
        let epDay = this.visitNamesAndDays.find(it => it.name === this.ep).day;
        let filteredDf = df.groupBy(df.columns.names())
        .where(`${LAB_VISIT_DAY} >= ${blDay} and ${LAB_VISIT_DAY} <= ${epDay} and ${LAB_TEST} = ${this.selectedLabValue}`)
        .aggregate();
        return filteredDf;
    }

}