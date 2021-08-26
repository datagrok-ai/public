import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { Filter, InputBase } from "datagrok-api/dg";
import $ from "cash-dom";
import { addColumnWithDrugPlusDosage } from "../data-preparation/data-preparation";
import { getUniqueValues } from "../data-preparation/utils";


export class PatientProfileView extends DG.ViewBase {

  options_lb_ae = {
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
          categories: [ "Basophils", "Urate", "Glucose" ], // fixed categories
          minField: 3,                    // min and max normal value 
          maxField: 4,                    // will be displayed with default color, otherwises "red"
          maxLimit: 5,                    // max number of categories
          alertColor: 'red',
        },
        markerShape: 'circle',
        height: '1flex',                  // height can be '30px', '20%', '3flex'
        show: 1,
      },

      // multi linechart 
      {
        tableName: 'patient_lb',
        type: 'line',
        multi: true,
        x: 'LBDY',
        y: 'LBSTRESN',
        splitByColumnName: 'LBTEST',                    // get categories from this column
        categories: [ "Basophils", "Urate", "Glucose" ],  // fixed categories
        maxLimit: 1,                                    // max number of linecharts 
        yType: 'value',
        markerShape: 'square',
        height: '1flex',
        show: 1,
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
      },

    ]

  }

  options_ex_cm = {
    series: [
      {
        tableName: 'patient_ex',
        title: 'Drug Exposure',
        type: 'timeLine',
        y: 'EXTRT_WITH_DOSE',
        x: [ 'EXSTDY', 'EXENDY' ],
        yType: 'value',
        color: 'red',
        markerShape: 'circle',
        height: '2flex',
        show: 1,
      },

      {
        tableName: 'patient_cm',
        title: 'Concomitant medication',
        type: 'timeLine',
        y: 'CMTRT',
        x: [ 'CMSTDY', 'CMENDY' ],
        yType: 'value',
        color: 'red',
        markerShape: 'circle',
        height: '2flex',
        show: 1,
      },

    ]

  }


  tableNames = [ 'lb', 'ae', 'ex', 'cm' ];
  tables = {};
  multiplot_lb_ae: any;
  multiplot_ex_cm: any;

  constructor() {
    super();

    let patientIds = Array.from(getUniqueValues(study.domains.dm, 'USUBJID'));
    let patienIdBoxPlot = ui.choiceInput('Patient Id', patientIds[ 0 ], patientIds);
    patienIdBoxPlot.onChanged((v) => {
      this.createTablesToAttach(patienIdBoxPlot.value);
      this.attachTablesToMultiplot(this.multiplot_lb_ae, this.options_lb_ae, [ 'lb', 'ae' ]);
      this.attachTablesToMultiplot(this.multiplot_ex_cm, this.options_ex_cm, [ 'ex', 'cm' ]);
    });

    this.createTablesToAttach(patientIds[ 0 ]);

    this.tables[ 'ae' ].plot.fromType('MultiPlot', {
      paramOptions: JSON.stringify(this.options_lb_ae),
    }).then((v: any) => {

      this.multiplot_lb_ae = v;
      this.attachTablesToMultiplot(this.multiplot_lb_ae, this.options_lb_ae, [ 'lb', 'ae' ]);
      this.tables[ 'ex' ].plot.fromType('MultiPlot', {
        paramOptions: JSON.stringify(this.options_ex_cm),
      }).then((v1: any) => {

        this.multiplot_ex_cm = v1;
        this.attachTablesToMultiplot(this.multiplot_ex_cm, this.options_ex_cm, [ 'ex', 'cm' ]);
        this.root.appendChild(
          ui.divV([ patienIdBoxPlot.root,
          ui.splitH([
            this.multiplot_lb_ae.root,
            this.multiplot_ex_cm.root
          ], { style: { width: '100%', height: '100%' } })
          ], { style: { width: '100%', height: '100%' } })

        );

      });

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

  private createTablesToAttach(myId: any) {
    this.tableNames.forEach(name => {
      this.tables[ name ] = study.domains[ name ].clone();
      this.tables[ name ].name = `patient_${name}`;
      this.tables[ name ].filter.init((i) => {
        let row = this.tables[ name ].row(i);
        return row[ 'USUBJID' ] === myId;
      })
    })
    this.tables[ 'ex' ] = addColumnWithDrugPlusDosage(this.tables[ 'ex' ], 'EXTRT', 'EXDOSE', 'EXDOSU', 'EXTRT_WITH_DOSE');
  }

}
