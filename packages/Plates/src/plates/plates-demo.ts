/* eslint-disable max-len */
import {PlateProperty, plateTemplates, plateTypes, savePlate, PlateTemplate} from './plates-crud';
import {
  initPlates, createAnalysisRun, saveAnalysisRunParameter,
  getOrCreateProperty, saveAnalysisResult,
} from './plates-crud';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {createPlateTemplate, createNewPlateForTemplate} from './plates-crud';
import {Plate} from '../plate/plate';
import {DrcAnalysis} from '../plate/analyses/drc/drc-analysis';
import {getDoseResponseSeries} from '../plate/analyses/drc/utils';
import {FIT_FUNCTION_4PL_REGRESSION, IFitSeries} from '@datagrok-libraries/statistics/src/fit/fit-curve';


export async function __createDummyPlateData() {
  await initPlates();
  await createDummyPlates();
  await createDummyPlatesFromExcel();

  await initPlates(true);
}

async function createDummyPlates() {
  let cellCountingTemplate = plateTemplates.find((t) => t.name === 'Cell counting');
  let createCcPlates = false;

  if (!cellCountingTemplate) {
    console.log('Creating "Cell counting" template...');
    cellCountingTemplate = await createPlateTemplate({
      name: 'Cell counting',
      description: 'Microscopy-based cell counting',
      plateProperties: [
        {name: 'Imaging device', choices: ['Kodak', 'Nikon'], type: DG.COLUMN_TYPE.STRING, is_required: true},
        {name: 'Status', choices: ['Pending', 'Filling', 'Measuring', 'Done'], type: DG.COLUMN_TYPE.STRING, is_required: true, default_value: 'Pending'},
        {name: 'Plate cell count', min: 0, max: 10000, type: DG.COLUMN_TYPE.INT, is_required: false}
      ],
      wellProperties: [
        {name: 'Well cell count', min: 0, max: 100, type: DG.COLUMN_TYPE.INT, is_required: false},
        {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], type: DG.COLUMN_TYPE.STRING, is_required: true}
      ]
    });
    createCcPlates = true;
  } else {
    console.log('"Cell counting" template found.');
    const ccPlatesCheck = await grok.data.db.query('Plates:Plts', `
      SELECT p.id FROM plts.plates p
      JOIN plts.plate_details pd ON p.id = pd.plate_id
      JOIN plts.properties prop ON pd.property_id = prop.id
      WHERE prop.name = 'Imaging device'
      LIMIT 1
    `);
    if (ccPlatesCheck.rowCount === 0) {
      console.log('No plates found for "Cell counting", creating them...');
      createCcPlates = true;
    }
  }

  let doseResponseTemplate = plateTemplates.find((t) => t.name === 'Dose-response');
  let createDrPlates = false;

  if (!doseResponseTemplate) {
    console.log('Creating "Dose-response" template...');
    doseResponseTemplate = await createPlateTemplate({
      name: 'Dose-response',
      description: 'Dose-response campaign',
      plateProperties: [
        {name: 'Project', type: DG.COLUMN_TYPE.STRING, choices: ['Modulators of mGluR5', 'Agonists for GPCR GPR139', 'Glutaminase Inhibitors for TNBC'], is_required: true},
        {name: 'Stage', type: DG.COLUMN_TYPE.STRING, choices: ['Lead generation', 'Lead optimization'], is_required: false},
        {name: 'Chemist', type: DG.COLUMN_TYPE.STRING, choices: ['John Marlowski', 'Mary Hopton'], is_required: true},
        {name: 'Biologist', type: DG.COLUMN_TYPE.STRING, choices: ['Anna Fei', 'Joan Dvorak'], is_required: true},
        {name: 'QC Passed', type: DG.COLUMN_TYPE.BOOL, is_required: true, default_value: 'false'},
        {name: 'Z-Score', min: 0, max: 3, type: DG.COLUMN_TYPE.FLOAT, is_required: false},
      ],
      wellProperties: [
        {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], type: DG.COLUMN_TYPE.STRING, is_required: false},
        {name: 'Role', choices: ['Control', 'Treatment'], type: DG.COLUMN_TYPE.STRING, is_required: true, default_value: 'Treatment'},
        {name: 'Concentration', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT, is_required: true},
        {name: 'Volume', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT, is_required: false},
        {name: 'Activity', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT, is_required: true},
      ]
    });
    createDrPlates = true;
  } else {
    console.log('"Dose-response" template found.');
    const drPlatesCheck = await grok.data.db.query('Plates:Plts', `
      SELECT p.id FROM plts.plates p
      JOIN plts.plate_details pd ON p.id = pd.plate_id
      JOIN plts.properties prop ON pd.property_id = prop.id
      WHERE prop.name = 'Project'
      LIMIT 1
    `);
    if (drPlatesCheck.rowCount === 0) {
      console.log('No plates found for "Dose-response", creating them...');
      createDrPlates = true;
    }
  }

  if (createCcPlates || createDrPlates) {
    await initPlates(true);
    cellCountingTemplate = plateTemplates.find((t) => t.name === 'Cell counting')!;
    doseResponseTemplate = plateTemplates.find((t) => t.name === 'Dose-response')!;
  }

  const getDemoValue = (property: PlateProperty) =>
    property.choices ? DG.Utils.random(typeof property.choices == 'string' ? JSON.parse(property.choices) : property.choices) :
      property.type === DG.COLUMN_TYPE.BOOL ? Math.random() > 0.5 :
        property.min !== undefined && property.max !== undefined ?
          property.type === DG.COLUMN_TYPE.INT ?
            Math.floor(property.min + Math.random() * (property.max - property.min)) :
            property.min + Math.random() * (property.max - property.min) :
          null;

  const templatesToProcess: [PlateTemplate | undefined, boolean][] = [
    [cellCountingTemplate, createCcPlates],
    [doseResponseTemplate, createDrPlates]
  ];

  for (const [template, shouldCreate] of templatesToProcess) {
    if (template && shouldCreate) {
      console.log(`Creating 20 demo plates for "${template.name}"...`);
      for (let i = 0; i < 20; i++) {
        const plate = await createNewPlateForTemplate(plateTypes[0], template);
        plate.details = Object.fromEntries(
          template.plateProperties.map((p) => [p.name!, getDemoValue(p as PlateProperty)])
        );

        for (const property of template.wellProperties)
          plate.data.col(property!.name!)?.init((_) => getDemoValue(property as PlateProperty));

        await savePlate(plate);

        if (template.name === 'Dose-response' && plate.id) {
          try {
            const drcAnalysis = new DrcAnalysis();
            const mappings = new Map<string, string>([
              ['Activity', 'Activity'],
              ['Concentration', 'Concentration'],
              ['SampleID', 'Sample'],
            ]);
            const params: Record<string, any> = {'Normalize': true};

            const sampleColName = mappings.get('SampleID')!;
            const concentrationColName = mappings.get('Concentration')!;
            const activityColName = mappings.get('Activity')!;

            const finalValueCol = activityColName;
            const controlColumns = ['High Control', 'Low Control'];
            if (params['Normalize'] && controlColumns.every((c) => plate.data.col(sampleColName)!.categories.includes(c))) {
              // ... normalization logic would go here, but forget it for the demo
            }

            const series = getDoseResponseSeries(plate, {
              value: finalValueCol,
              concentration: concentrationColName,
              groupBy: sampleColName,
            });

            const seriesVals = Object.entries(series);
            if (seriesVals.length === 0 || !seriesVals.some(([_, s]) => s.points.length > 1)) {
              console.warn(`DRC: No series data for demo plate ${plate.id}, skipping analysis save.`);
              continue;
            }

            const roleCol = DG.Column.string(sampleColName, seriesVals.length).init((i) => seriesVals[i][0]);
            const curveCol = DG.Column.string('Curve', seriesVals.length);
            curveCol.semType = 'fit';

            curveCol.init((i) => {
              const currentSeries = seriesVals[i][1];
              const seriesData: IFitSeries = {
                name: currentSeries.name,
                points: currentSeries.points,
                fitFunction: FIT_FUNCTION_4PL_REGRESSION,
                parameters: undefined,
                clickToToggle: true,
                droplines: ['IC50'],
                showPoints: 'points',
              };
              const chartData = {
                chartOptions: {
                  xAxisName: concentrationColName,
                  yAxisName: finalValueCol,
                  logX: true,
                  title: `${seriesVals[i][0]}`,
                },
                series: [seriesData],
              };
              return JSON.stringify(chartData);
            });

            const resultsDf = DG.DataFrame.fromColumns([roleCol, curveCol]);

            const statsToAdd: Record<string, string> = {
              'interceptX': 'IC50', 'slope': 'Hill Slope', 'rSquared': 'R Squared',
              'bottom': 'Min', 'top': 'Max', 'auc': 'AUC',
            };
            const outputNames = new Set(drcAnalysis.outputs.map((o) => o.name));

            for (const [statName, colName] of Object.entries(statsToAdd)) {
              if (outputNames.has(colName)) {
                const col = await grok.functions.call('Curves:AddStatisticsColumn', {
                  table: resultsDf,
                  colName: 'Curve',
                  propName: statName,
                  seriesNumber: 0,
                });
                col.name = colName;
              }
            }

            const groupColumn = resultsDf.columns.byIndex(0).name;
            const groups = resultsDf.col(groupColumn)!.categories;

            const runId = await createAnalysisRun(plate.id, drcAnalysis.name, groups);

            for (const param of drcAnalysis.parameters) {
              await saveAnalysisRunParameter({
                runId: runId,
                propertyName: param.name,
                propertyType: param.type as DG.TYPE,
                value: params[param.name] ?? param.defaultValue,
              });
            }

            const outputProperties = new Map<string, PlateProperty>();
            for (const output of drcAnalysis.outputs) {
              const prop = await getOrCreateProperty(output.name, output.type as DG.TYPE, 'plate');
              outputProperties.set(output.name, prop);
            }

            for (let rowIdx = 0; rowIdx < resultsDf.rowCount; rowIdx++) {
              const row = resultsDf.row(rowIdx);
              const groupName = row.get(groupColumn);
              const groupCombination = [groupName];

              for (const output of drcAnalysis.outputs) {
                if (resultsDf.columns.contains(output.name)) {
                  const value = row.get(output.name);
                  const prop = outputProperties.get(output.name)!;

                  if (value !== null && value !== undefined) {
                    await saveAnalysisResult({
                      runId: runId,
                      propertyId: prop.id,
                      propertyName: prop.name,
                      propertyType: prop.type,
                      value: value,
                      groupCombination: groupCombination,
                    });
                  }
                }
              }
            }
            console.log(`Saved demo DRC analysis (runId: ${runId}) for plate ${plate.id}`);
          } catch (e) {
            console.error(`Failed to create demo analysis for plate ${plate.id}:`, e);
          }
        }
      }
    } else if (template) {
      console.log(`Plates for "${template.name}" already exist. Skipping.`);
    }
  }
}

async function createDummyPlatesFromExcel() {
  let excelPlateTemplate = plateTemplates.find((t) => t.name === 'Excel');
  let createExcelPlates = false;

  if (!excelPlateTemplate) {
    console.log('Creating "Excel" template...');
    excelPlateTemplate = await createPlateTemplate({
      name: 'Excel',
      description: 'Excel plates',
      plateProperties: [{name: 'Barcode', type: DG.COLUMN_TYPE.STRING}],
      wellProperties: [
        {name: 'raw data', type: DG.COLUMN_TYPE.FLOAT},
        {name: 'plate layout', type: DG.COLUMN_TYPE.STRING},
        {name: 'concentrations', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT}
      ]
    });
    createExcelPlates = true;
  } else {
    console.log('"Excel" template found.');
    const excelPlatesCheck = await grok.data.db.query('Plates:Plts', `
      SELECT p.id FROM plts.plates p
      JOIN plts.plate_well_values pwv ON p.id = pwv.plate_id
      JOIN plts.properties prop ON pwv.property_id = prop.id
      WHERE prop.name = 'raw data'
      LIMIT 1
    `);
    if (excelPlatesCheck.rowCount === 0) {
      console.log('No plates found for "Excel", creating them...');
      createExcelPlates = true;
    }
  }

  if (createExcelPlates) {
    await initPlates(true);
    excelPlateTemplate = plateTemplates.find((t) => t.name === 'Excel')!;

    console.log('Creating demo plates from Excel files...');
    const xlsxFiles = await grok.dapi.files.list('System.DemoFiles/hts/xlsx_plates', false, '*.xlsx');
    for (const file of xlsxFiles) {
      const plate = await Plate.fromExcelFile(file);
      plate.plateTemplateId = excelPlateTemplate.id;
      plate.details = {Barcode: (Math.round(Math.random() * 100000)).toString().padStart(8, '0')};
      await savePlate(plate, {autoCreateProperties: false});
    }
  } else {
    console.log('Plates for "Excel" template already exist. Skipping.');
  }
}
