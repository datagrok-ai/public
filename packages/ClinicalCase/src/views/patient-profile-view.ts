import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { AE_END_DAY, AE_START_DAY, AE_TERM, CON_MED_END_DAY, CON_MED_NAME, CON_MED_START_DAY, INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, INV_DRUG_END_DAY, INV_DRUG_NAME, INV_DRUG_START_DAY, LAB_DAY, LAB_HI_LIM_N, LAB_LO_LIM_N, LAB_RES_N, LAB_TEST, SUBJECT_ID } from "../columns-constants";
import { requiredColumnsByView } from "../constants";
import { addColumnWithDrugPlusDosage, dynamicComparedToBaseline, labDynamicComparedToMinMax, labDynamicRelatedToRef } from "../data-preparation/data-preparation";
import { getUniqueValues } from "../data-preparation/utils";
import { ILazyLoading } from "../lazy-loading/lazy-loading";
import { _package } from "../package";
import { checkMissingDomains } from "./utils";


export class PatientProfileView extends DG.ViewBase implements ILazyLoading {

  options_lb_ae_ex_cm = {
    series: [
      {
        tableName: 'patient_lb',
        title: 'Lab values',
        type: 'scatter',
        x: LAB_DAY,
        y: LAB_TEST,
        // extraFields is an array to load into echart data arrays
        // all fields later combined into one array [x, y, extraFields]
        // user can address fields by index, for instance index 3 means field "LBORNRLO"
        extraFields: [LAB_RES_N, LAB_LO_LIM_N, LAB_HI_LIM_N],
        yType: 'category',                // can be 'value' or 'category'
        statusChart: {
          valueField: 2,                  // index of field with test value
          splitByColumnName: LAB_RES_N,    // column to get categories
          categories: [''], // fixed categories
          minField: 3,                    // min and max normal value 
          maxField: 4,                    // will be displayed with default color, otherwises "red"
          alertColor: 'red',
        },
        markerShape: 'circle',
        height: '1flex',                  // height can be '30px', '20%', '3flex'
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate',
        multiEdit: {
          options: study.domains.lb ? Array.from(getUniqueValues(study.domains.lb, LAB_TEST)) : [],
          selectedValues: [],
          editValue: 'category',
          updateTitle: false
        }
      },

      //linechart 
      {
        tableName: 'patient_lb',
        title: 'Lab values line chart',
        type: 'line',
        multiLineFieldIndex: 2, //index of field by which to split multiple graphs
        x: LAB_DAY,
        y: 'LAB_DYNAMIC_BL',
        extraFields: [LAB_TEST, LAB_RES_N, LAB_LO_LIM_N, LAB_HI_LIM_N, `BL_${LAB_RES_N}`, `min(${LAB_RES_N})`, `max(${LAB_RES_N})`],
        splitByColumnName: LAB_TEST,                    // get categories from this column
        categories: [''],  // fixed categories
        maxLimit: 1,                                    // max number of linecharts 
        yType: 'category',
        markerShape: 'square',
        height: '1flex',
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate',
        multiEdit: {
          options: study.domains.lb ? Array.from(getUniqueValues(study.domains.lb, LAB_TEST)) : [],
          selectedValues: [],
          editValue: 'category',
          updateTitle: false
        },
        comboEdit: {
          options: ['From BL', 'Min/Max', 'Related to ref'],
          selectedValue: 'From BL',
          editValue: 'y',
          editName: 'Type',
          values: { 'From BL': 'LAB_DYNAMIC_BL', 'Min/Max': 'LAB_DYNAMIC_MIN_MAX', 'Related to ref': 'LAB_DYNAMIC_REF' },
          additionalParams: {
            'Related to ref': {
              markLine: {
                lineStyle: {
                  color: '#333'
                },
                data: [
                  {
                    name: 'lower lim',
                    yAxis: -1,
                    label: {
                      formatter: '{b}',
                      position: 'start'
                    }
                  },
                  {
                    name: 'upper lim',
                    yAxis: 1,
                    label: {
                      formatter: '{b}',
                      position: 'start'
                    }
                  },
                ]
              },
              min: -1.5,
              max: 1.5
            }
          }
        },
        showLegend: true,
      },

      // timeLines
      {
        tableName: 'patient_ae',
        title: 'Adverse Events',
        type: 'timeLine',
        y: AE_TERM,                      // category column
        x: [AE_START_DAY, AE_END_DAY],          // [startTime, endTime]
        yType: 'category',
        color: 'red',                     // color of marker
        markerShape: 'circle',
        height: '2flex',
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate'
      },

      {
        tableName: 'patient_ex',
        title: 'Drug Exposure',
        type: 'timeLine',
        y: 'EXTRT_WITH_DOSE',
        x: [INV_DRUG_START_DAY, INV_DRUG_END_DAY],
        yType: 'category',
        color: 'red',
        markerShape: 'circle',
        height: '2flex',
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate'
      },

      {
        tableName: 'patient_cm',
        title: 'Concomitant medication',
        type: 'timeLine',
        y: CON_MED_NAME,
        x: [CON_MED_START_DAY, CON_MED_END_DAY],
        yType: 'category',
        color: 'red',
        markerShape: 'circle',
        height: '2flex',
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate'
      },

    ]

  }

  tableNamesAndFields = {
    'lb': { 'start': LAB_DAY },
    'ae': { 'start': AE_START_DAY, 'end': AE_END_DAY },
    'ex': { 'start': INV_DRUG_START_DAY, 'end': INV_DRUG_END_DAY },
    'cm': { 'start': CON_MED_START_DAY, 'end': CON_MED_END_DAY }
  };
  tables = {};
  multiplot_lb_ae_ex_cm: any;

  constructor(name) {
    super({});
    this.name = name;
    this.helpUrl = `${_package.webRoot}/views_help/patient_profile.md`;
  }

  loaded: boolean;

  load(): void {
    checkMissingDomains(requiredColumnsByView[this.name], this);
  }

  createView(): void {
    let patientIds = Array.from(getUniqueValues(study.domains.dm, SUBJECT_ID));
    let patienIdBoxPlot = ui.choiceInput('', patientIds[0], patientIds);
    patienIdBoxPlot.onChanged((v) => {
      this.updateTablesToAttach(patienIdBoxPlot.value);
      this.attachTablesToMultiplot(this.multiplot_lb_ae_ex_cm, this.options_lb_ae_ex_cm, ['lb', 'ae', 'ex', 'cm']);
      this.multiplot_lb_ae_ex_cm.setAdditionalParams(1, this.options_lb_ae_ex_cm.series[1].comboEdit.additionalParams[this.options_lb_ae_ex_cm.series[1].comboEdit.selectedValue])
      let updateTitle = this.options_lb_ae_ex_cm.series[1].multiEdit.selectedValues.length ? true : false;
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(1, this.options_lb_ae_ex_cm.series[1].multiEdit.selectedValues, updateTitle);
      this.multiplot_lb_ae_ex_cm.updatePlotByYAxis(1, this.options_lb_ae_ex_cm.series[1].comboEdit.values[this.options_lb_ae_ex_cm.series[1].comboEdit.selectedValue]);
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(0, this.options_lb_ae_ex_cm.series[0].multiEdit.selectedValues, false);
    });

    this.updateTablesToAttach(patientIds[0]);


    this.tables['ae'].plot.fromType('MultiPlot', {
      paramOptions: JSON.stringify(this.options_lb_ae_ex_cm),
    }).then((v: any) => {

      this.multiplot_lb_ae_ex_cm = v;
      this.attachTablesToMultiplot(this.multiplot_lb_ae_ex_cm, this.options_lb_ae_ex_cm, ['lb', 'ae', 'ex', 'cm']);
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(0, this.options_lb_ae_ex_cm.series[0].multiEdit.selectedValues, false); //to clear scattr plot after creation
      this.setRibbonPanels([
        [
          ui.iconFA('chevron-left', () => {
            //@ts-ignore
            let current = patienIdBoxPlot.input.selectedIndex;
            if (current != 0) {
              current--;
              //@ts-ignore
              patienIdBoxPlot.value = patienIdBoxPlot.input.options[current].value;
              //@ts-ignore
              patienIdBoxPlot.input.selectedIndex = current
              patienIdBoxPlot.fireChanged();
            }
          }),
        ],
        [
          patienIdBoxPlot.root
        ],
        [
          ui.iconFA('chevron-right', () => {
            //@ts-ignore
            let current = patienIdBoxPlot.input.selectedIndex;
            //@ts-ignore
            let length = patienIdBoxPlot.input.length;
            if (current != length) {
              current++;
              //@ts-ignore
              patienIdBoxPlot.value = patienIdBoxPlot.input.options[current].value;
              //@ts-ignore
              patienIdBoxPlot.input.selectedIndex = current
              patienIdBoxPlot.fireChanged();
            }
          })
        ]
      ]);
      this.root.className = 'grok-view ui-box';
      this.root.appendChild(this.multiplot_lb_ae_ex_cm.root);
    });
  }

  private attachTablesToMultiplot(plot: any, options: any, tableNames: string[]) {
    const tablesToAttach = {};
    tableNames.forEach(name => {
      tablesToAttach[`patient_${name}`] = this.tables[name]
    });
    plot.tables = tablesToAttach;
    plot.options = options;
    plot.onTableAttached();
  }


  /*   private createTablesToAttach() {
      Object.keys(this.tableNamesAndFields).forEach(name => {
        this.tables[ name ] = study.domains[ name ].clone();
        this.tables[ name ].name = `patient_${name}`;
      })
    } */

  private updateTablesToAttach(myId: any) {
    Object.keys(this.tableNamesAndFields).forEach(name => {
      this.tables[name] = study.domains[name].clone()
        .groupBy(study.domains[name].columns.names())
        .where(`${SUBJECT_ID} = ${myId}`)
        .aggregate();
      this.tables[name].name = `patient_${name}`;
    });
    this.options_lb_ae_ex_cm['xAxisMinMax'] = this.extractMinAndMaxValuesForXAxis();
    addColumnWithDrugPlusDosage(this.tables['ex'], INV_DRUG_NAME, INV_DRUG_DOSE, INV_DRUG_DOSE_UNITS, 'EXTRT_WITH_DOSE');
    dynamicComparedToBaseline(this.tables['lb'], LAB_TEST, LAB_RES_N, this.options_lb_ae_ex_cm['xAxisMinMax']['minX'], LAB_DAY, 'LAB_DYNAMIC_BL', false);
    labDynamicComparedToMinMax(this.tables['lb'], 'LAB_DYNAMIC_MIN_MAX');
    labDynamicRelatedToRef(this.tables['lb'], 'LAB_DYNAMIC_REF');
  }

  private extractMinAndMaxValuesForXAxis() {
    let min = null;
    let max = null;
    Object.keys(this.tableNamesAndFields).forEach(tableName => {
      let table = this.tables[tableName];
      let minColName = this.tableNamesAndFields[tableName]['start'];
      let maxColName = this.tableNamesAndFields[tableName]['end'] ?? this.tableNamesAndFields[tableName]['start'];
      let newMin = table.getCol(minColName).stats['min'];
      min = min !== null || newMin > min ? min : newMin;
      let newMax = table.getCol(maxColName).stats['max'];
      max = max !== null || newMax < max ? max : newMax;
    })
    return { minX: min, maxX: max };
  }

}
