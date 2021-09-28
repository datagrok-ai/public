import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { addColumnWithDrugPlusDosage, labDynamicComparedToBaseline, labDynamicComparedToMinMax } from "../data-preparation/data-preparation";
import { getUniqueValues } from "../data-preparation/utils";


export class PatientProfileView extends DG.ViewBase {

  options_lb_ae_ex_cm = {
    series: [
      {
        tableName: 'patient_lb',
        title: 'Lab values',
        type: 'scatter',
        x: 'LBDY',
        y: 'LBTEST',
        // extraFields is an array to load into echart data arrays
        // all fields later combined into one array [x, y, extraFields]
        // user can address fields by index, for instance index 3 means field "LBORNRLO"
        extraFields: [ 'LBORRES', 'LBORNRLO', 'LBORNRHI' ],
        yType: 'category',                // can be 'value' or 'category'
        statusChart: {
          valueField: 2,                  // index of field with test value
          splitByColumnName: 'LBTEST',    // column to get categories
          categories: [ '' ], // fixed categories
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
          options: Array.from(getUniqueValues(study.domains.lb, 'LBTEST')), 
          selectedValues: [], 
          editValue: 'category'}
      },

       //linechart 
      {
        tableName: 'patient_lb',
        title: 'Lab values line chart',
        type: 'line',
        multiLineFieldIndex: 2, //index of field by which to split multiple graphs
        x: 'LBDY',
        y: 'LAB_DYNAMIC_BL',
        extraFields: [ 'LBTEST', 'LBORRES', 'LBORNRLO', 'LBORNRHI', 'BL_LBORRES', 'min(LBORRES_FLOAT)', 'max(LBORRES_FLOAT)' ],
        splitByColumnName: 'LBTEST',                    // get categories from this column
        categories: [ '' ],  // fixed categories
        maxLimit: 1,                                    // max number of linecharts 
        yType: 'category',
        markerShape: 'square',
        height: '1flex',
        show: 1,
        yLabelWidth: 50,
        yLabelOverflow: 'truncate',
        multiEdit: {
          options: Array.from(getUniqueValues(study.domains.lb, 'LBTEST')), 
          selectedValues: [], 
          editValue: 'category'},
        comboEdit: {
          options: ['From BL', 'Min/Max'], 
          selectedValue: 'From BL', 
          editValue: 'y', 
          editName: 'Type', 
          values: {'From BL': 'LAB_DYNAMIC_BL', 'Min/Max': 'LAB_DYNAMIC_MIN_MAX'}},
        showLegend: true
      }, 

      // timeLines
      {
        tableName: 'patient_ae',
        title: 'Adverse Events',
        type: 'timeLine',
        y: 'AETERM',                      // category column
        x: [ 'AESTDY', 'AEENDY' ],          // [startTime, endTime]
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
        x: [ 'EXSTDY', 'EXENDY' ],
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
        y: 'CMTRT',
        x: [ 'CMSTDY', 'CMENDY' ],
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

  tableNamesAndFields =  {
    'lb': {'start': 'LBDY'}, 
    'ae': {'start': 'AESTDY', 'end': 'AEENDY'}, 
    'ex': {'start': 'EXSTDY', 'end': 'EXENDY'},
    'cm': {'start': 'CMSTDY', 'end': 'CMENDY'}
  };
  tables = {};
  multiplot_lb_ae_ex_cm: any;

  constructor(name) {
    super(name);

    this.name = name;

    let patientIds = Array.from(getUniqueValues(study.domains.dm, 'USUBJID'));
    let patienIdBoxPlot = ui.choiceInput('', patientIds[ 0 ], patientIds);
    patienIdBoxPlot.onChanged((v) => {
      this.updateTablesToAttach(patienIdBoxPlot.value);
      this.attachTablesToMultiplot(this.multiplot_lb_ae_ex_cm, this.options_lb_ae_ex_cm, [ 'lb', 'ae', 'ex', 'cm' ]);
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(1, this.options_lb_ae_ex_cm.series[1].multiEdit.selectedValues, true);
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(0, this.options_lb_ae_ex_cm.series[0].multiEdit.selectedValues, false);
    });


    this.createTablesToAttach();
    this.updateTablesToAttach(patientIds[ 0 ]);


    this.tables[ 'ae' ].plot.fromType('MultiPlot', {
      paramOptions: JSON.stringify(this.options_lb_ae_ex_cm),
    }).then((v: any) => {

      this.multiplot_lb_ae_ex_cm = v;
      this.attachTablesToMultiplot(this.multiplot_lb_ae_ex_cm, this.options_lb_ae_ex_cm, [ 'lb', 'ae', 'ex', 'cm' ]);
      this.multiplot_lb_ae_ex_cm.updatePlotByCategory(0, this.options_lb_ae_ex_cm.series[0].multiEdit.selectedValues, false); //to clear scattr plot after creation
      this.setRibbonPanels([
        [
          ui.iconFA('chevron-left',()=>{
            //@ts-ignore
            let current = patienIdBoxPlot.input.selectedIndex;
            if (current!=0){
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
          ui.iconFA('chevron-right',()=>{
            //@ts-ignore
            let current = patienIdBoxPlot.input.selectedIndex;
            //@ts-ignore
            let length = patienIdBoxPlot.input.length;
            if (current!=length){
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
      tablesToAttach[ `patient_${name}` ] = this.tables[ name ]
    });
    plot.tables = tablesToAttach;
    plot.options = options;
    plot.onTableAttached();
  }


  private createTablesToAttach() {
    Object.keys(this.tableNamesAndFields).forEach(name => {
      this.tables[ name ] = study.domains[ name ].clone();
      this.tables[ name ].name = `patient_${name}`;
    })
  }

  private updateTablesToAttach(myId: any){
    Object.keys(this.tableNamesAndFields).forEach(name => {
      this.tables[ name ].filter.init((i) => {
        let row = this.tables[ name ].row(i);
        return row[ 'USUBJID' ] === myId;
      })
    })
    addColumnWithDrugPlusDosage(this.tables[ 'ex' ], 'EXTRT', 'EXDOSE', 'EXDOSU', 'EXTRT_WITH_DOSE');
    this.options_lb_ae_ex_cm['xAxisMinMax'] = this.extractMinAndMaxValuesForXAxis();
    labDynamicComparedToBaseline(this.tables[ 'lb' ],  this.options_lb_ae_ex_cm['xAxisMinMax']['minX'], 'LAB_DYNAMIC_BL');
    labDynamicComparedToMinMax(this.tables[ 'lb' ], 'LAB_DYNAMIC_MIN_MAX');
  }

  private extractMinAndMaxValuesForXAxis() {
    let min = null;
    let max = null;
    Object.keys(this.tables).forEach(table => {
      let minColName = this.tableNamesAndFields[ table ][ 'start' ];
      let maxColName = this.tableNamesAndFields[ table ][ 'end' ] ?? this.tableNamesAndFields[ table ][ 'start' ];
      let newMin = this.tables[ table ].getCol(minColName).stats[ 'min' ];
      min = min !== null || newMin > min ? min : newMin;
      let newMax = this.tables[ table ].getCol(maxColName).stats[ 'max' ];
      max = max !== null || newMax < max ? max : newMax;
    })
    return {minX: min, maxX: max};
  }

}
