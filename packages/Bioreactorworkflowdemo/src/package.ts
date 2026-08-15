/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {IRuntimeLinkController, PipelineConfiguration} from '@datagrok-libraries/compute-api';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import timezone from 'dayjs/plugin/timezone';
import {simulateDay} from './bioreactor-model';
export * from './package.g';

dayjs.extend(utc);
dayjs.extend(timezone);

export const _package = new DG.Package();

const DAILY_SUMMARY_COLUMNS = [
  {field: 'solutionAdded', name: 'Solution added', units: 'L'},
  {field: 'substrateAdded', name: 'Substrate added', units: 'g'},
  {field: 'incomingVolume', name: 'Incoming volume', units: 'L'},
  {field: 'incomingBiomass', name: 'Incoming biomass', units: 'g/L'},
  {field: 'incomingSubstrate', name: 'Incoming substrate', units: 'g/L'},
  {field: 'maximumGrowthRate', name: 'Maximum growth rate', units: '1/h'},
  {field: 'halfSaturationConstant', name: 'Monod constant', units: 'g/L'},
  {field: 'yieldCoefficient', name: 'Biomass yield', units: 'g/g'},
  {field: 'decayConstant', name: 'Decay rate', units: '1/h'},
  {field: 'dayDuration', name: 'Day duration', units: 'h'},
  {field: 'finalVolume', name: 'Final volume', units: 'L'},
  {field: 'finalBiomass', name: 'Final biomass', units: 'g/L'},
  {field: 'finalSubstrate', name: 'Final substrate', units: 'g/L'},
  {field: 'finalBiomassMass', name: 'Final biomass mass', units: 'g'},
  {field: 'finalSubstrateMass', name: 'Final substrate mass', units: 'g'},
  {field: 'biomassMassChange', name: 'Biomass mass change', units: 'g'},
  {field: 'substrateConsumed', name: 'Substrate consumed', units: 'g'},
] as const;

function makeDailySummary(dayGetter: (field: string) => unknown[]): DG.DataFrame {
  const firstField = DAILY_SUMMARY_COLUMNS[0].field;
  const dayCount = dayGetter(firstField).length;
  const columns: DG.Column[] = [
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Day', Array.from({length: dayCount}, (_, index) => index + 1)),
  ];

  for (const {field, name, units} of DAILY_SUMMARY_COLUMNS) {
    const column = DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, dayGetter(field));
    column.meta.units = units;
    columns.push(column);
  }

  const summary = DG.DataFrame.fromColumns(columns);
  summary.name = 'Cultivation summary';
  return summary;
}

function populateDailySummary({controller}: {controller: IRuntimeLinkController}): void {
  const summary = makeDailySummary((field) => controller.getAll(`day_${field}`) ?? []);
  controller.setAll('summary', summary, 'restricted');
}

//name: BioreactorWorkflow
//description: Configure and simulate daily bacterial growth in a bioreactor
//meta.role: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function bioreactorWorkflow(params: {version?: string}): PipelineConfiguration {
  void params;

  return {
    id: 'bioreactorWorkflow',
    nqName: 'BioreactorWorkflowDemo:BioreactorWorkflow',
    version: '1.0',
    type: 'static',
    friendlyName: 'Bioreactor workflow demo',
    steps: [
      {
        id: 'bioreactorConfiguration',
        nqName: 'BioreactorWorkflowDemo:BioreactorConfiguration',
        friendlyName: 'Bioreactor configuration',
      },
      {
        id: 'dailyCultivation',
        type: 'dynamic',
        friendlyName: 'Daily cultivation',
        stepTypes: [{
          id: 'dayCalculation',
          nqName: 'BioreactorWorkflowDemo:DayCalculation',
          friendlyName: 'Day',
        }],
        initialSteps: [{id: 'dayCalculation'}],
        links: [
          {
            id: 'chain-volume',
            base: 'base:expand(dayCalculation)',
            from: 'previous:same(@base, dayCalculation)/finalVolume',
            to: 'next:after+(@base, dayCalculation)/incomingVolume',
            defaultRestrictions: {next: 'disabled'},
          },
          {
            id: 'chain-biomass',
            base: 'base:expand(dayCalculation)',
            from: 'previous:same(@base, dayCalculation)/finalBiomass',
            to: 'next:after+(@base, dayCalculation)/incomingBiomass',
            defaultRestrictions: {next: 'disabled'},
          },
          {
            id: 'chain-substrate',
            base: 'base:expand(dayCalculation)',
            from: 'previous:same(@base, dayCalculation)/finalSubstrate',
            to: 'next:after+(@base, dayCalculation)/incomingSubstrate',
            defaultRestrictions: {next: 'disabled'},
          },
        ],
      },
      {
        id: 'summary',
        nqName: 'BioreactorWorkflowDemo:Summary',
        friendlyName: 'Summary',
      },
    ],
    links: [
      {
        id: 'initial-volume',
        from: 'source:bioreactorConfiguration/volume',
        to: 'target:dailyCultivation/first(dayCalculation)/incomingVolume',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'initial-biomass',
        from: 'source:bioreactorConfiguration/biomass',
        to: 'target:dailyCultivation/first(dayCalculation)/incomingBiomass',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'initial-substrate',
        from: 'source:bioreactorConfiguration/substrate',
        to: 'target:dailyCultivation/first(dayCalculation)/incomingSubstrate',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'maximum-growth-rate',
        from: 'source:bioreactorConfiguration/maximumGrowthRate',
        to: 'target:dailyCultivation/all(dayCalculation)/maximumGrowthRate',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'half-saturation-constant',
        from: 'source:bioreactorConfiguration/halfSaturationConstant',
        to: 'target:dailyCultivation/all(dayCalculation)/halfSaturationConstant',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'yield-coefficient',
        from: 'source:bioreactorConfiguration/yieldCoefficient',
        to: 'target:dailyCultivation/all(dayCalculation)/yieldCoefficient',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'decay-constant',
        from: 'source:bioreactorConfiguration/decayConstant',
        to: 'target:dailyCultivation/all(dayCalculation)/decayConstant',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'day-duration',
        from: 'source:bioreactorConfiguration/dayDuration',
        to: 'target:dailyCultivation/all(dayCalculation)/dayDuration',
        defaultRestrictions: {target: 'disabled'},
      },
      {
        id: 'daily-summary',
        from: `day_(template):dailyCultivation/all(dayCalculation)/ ${
          DAILY_SUMMARY_COLUMNS.map(({field}) => field).join(' | ')
        }`,
        to: 'summary:summary/summaryInput',
        handler: populateDailySummary,
      },
    ],
  };
}

//name: BioreactorConfiguration
//description: Define the initial reactor state and bacterial growth parameters
//input: double initialVolume = 1 {caption: Initial volume; units: L; category: Initial state; min: 0.001} [Initial liquid volume in the reactor]
//input: double initialBiomass = 0.1 {caption: Initial biomass; units: g/L; category: Initial state; min: 0} [Initial biomass concentration]
//input: double initialSubstrate = 10 {caption: Initial substrate; units: g/L; category: Initial state; min: 0} [Initial limiting-substrate concentration]
//input: double muMax = 0.4 {caption: Maximum growth rate; units: 1/h; category: Growth kinetics; min: 0.001} [Maximum specific growth rate]
//input: double monodKs = 0.1 {caption: Monod constant; units: g/L; category: Growth kinetics; min: 0.001} [Substrate concentration at half of the maximum growth rate]
//input: double biomassYield = 0.5 {caption: Biomass yield; units: g/g; category: Growth kinetics; min: 0.001} [Biomass formed per unit of consumed substrate]
//input: double decayRate = 0.01 {caption: Decay rate; units: 1/h; category: Growth kinetics; min: 0} [First-order biomass decay rate]
//input: double maxVolume = 5 {caption: Maximum volume; units: L; category: Reactor; min: 0.001} [Reactor working-volume limit]
//output: double volume {caption: Initial volume; units: L; category: Initial state}
//output: double biomass {caption: Initial biomass; units: g/L; category: Initial state}
//output: double substrate {caption: Initial substrate; units: g/L; category: Initial state}
//output: double maximumGrowthRate {caption: Maximum growth rate; units: 1/h; category: Growth kinetics}
//output: double halfSaturationConstant {caption: Monod constant; units: g/L; category: Growth kinetics}
//output: double yieldCoefficient {caption: Biomass yield; units: g/g; category: Growth kinetics}
//output: double decayConstant {caption: Decay rate; units: 1/h; category: Growth kinetics}
//output: double dayDuration {caption: Day duration; units: h; category: Simulation}
//output: double workingVolumeLimit {caption: Maximum volume; units: L; category: Reactor}
export function bioreactorConfiguration(
  initialVolume: number,
  initialBiomass: number,
  initialSubstrate: number,
  muMax: number,
  monodKs: number,
  biomassYield: number,
  decayRate: number,
  maxVolume: number,
) {
  return {
    volume: initialVolume,
    biomass: initialBiomass,
    substrate: initialSubstrate,
    maximumGrowthRate: muMax,
    halfSaturationConstant: monodKs,
    yieldCoefficient: biomassYield,
    decayConstant: decayRate,
    dayDuration: 24,
    workingVolumeLimit: maxVolume,
  };
}

//name: DayCalculation
//description: Apply daily feeding and simulate bacterial growth for one cultivation day
//input: double solutionAdded = 0 {caption: Solution added; units: L; category: Daily feeding; min: 0} [Substrate-free solution added at the start of the day]
//input: double substrateAdded = 0 {caption: Substrate added; units: g; category: Daily feeding; min: 0} [Dry limiting substrate added at the start of the day]
//input: double incomingVolume = 1 {caption: Incoming volume; units: L; category: Incoming state; min: 0.001} [Final liquid volume from the preceding day]
//input: double incomingBiomass = 0.1 {caption: Incoming biomass; units: g/L; category: Incoming state; min: 0} [Final biomass concentration from the preceding day]
//input: double incomingSubstrate = 10 {caption: Incoming substrate; units: g/L; category: Incoming state; min: 0} [Final substrate concentration from the preceding day]
//input: double maximumGrowthRate = 0.4 {caption: Maximum growth rate; units: 1/h; category: Growth kinetics; min: 0.001}
//input: double halfSaturationConstant = 0.1 {caption: Monod constant; units: g/L; category: Growth kinetics; min: 0.001}
//input: double yieldCoefficient = 0.5 {caption: Biomass yield; units: g/g; category: Growth kinetics; min: 0.001}
//input: double decayConstant = 0.01 {caption: Decay rate; units: 1/h; category: Growth kinetics; min: 0}
//input: double dayDuration = 24 {caption: Day duration; units: h; category: Simulation; min: 0}
//output: double finalVolume {caption: Final volume; units: L; category: Final state}
//output: double finalBiomass {caption: Final biomass; units: g/L; category: Final state}
//output: double finalSubstrate {caption: Final substrate; units: g/L; category: Final state}
//output: double finalBiomassMass {caption: Final biomass mass; units: g; category: Mass balance}
//output: double finalSubstrateMass {caption: Final substrate mass; units: g; category: Mass balance}
//output: double biomassMassChange {caption: Biomass mass change; units: g; category: Mass balance}
//output: double substrateConsumed {caption: Substrate consumed; units: g; category: Mass balance}
//output: dataframe dailyProfile {caption: Daily profile}
export function dayCalculation(
  solutionAdded: number,
  substrateAdded: number,
  incomingVolume: number,
  incomingBiomass: number,
  incomingSubstrate: number,
  maximumGrowthRate: number,
  halfSaturationConstant: number,
  yieldCoefficient: number,
  decayConstant: number,
  dayDuration: number,
) {
  const result = simulateDay(
    solutionAdded,
    substrateAdded,
    incomingVolume,
    incomingBiomass,
    incomingSubstrate,
    maximumGrowthRate,
    halfSaturationConstant,
    yieldCoefficient,
    decayConstant,
    dayDuration,
  );
  const dailyProfile = DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Time', result.profile.map((point) => point.time)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Volume', result.profile.map((point) => point.volume)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Biomass', result.profile.map((point) => point.biomass)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Substrate', result.profile.map((point) => point.substrate)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Biomass mass', result.profile.map((point) => point.biomassMass)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'Substrate mass', result.profile.map((point) => point.substrateMass)),
  ]);
  const units: Record<string, string> = {
    Time: 'h',
    Volume: 'L',
    Biomass: 'g/L',
    Substrate: 'g/L',
    'Biomass mass': 'g',
    'Substrate mass': 'g',
  };
  for (const [columnName, unit] of Object.entries(units))
    dailyProfile.col(columnName)!.meta.units = unit;

  return {
    finalVolume: result.finalVolume,
    finalBiomass: result.finalBiomass,
    finalSubstrate: result.finalSubstrate,
    finalBiomassMass: result.finalBiomassMass,
    finalSubstrateMass: result.finalSubstrateMass,
    biomassMassChange: result.biomassMassChange,
    substrateConsumed: result.substrateConsumed,
    dailyProfile,
  };
}

//name: Summary
//description: Collect all daily cultivation inputs and outputs in one table
//input: dataframe summaryInput {caption: Daily cultivation data}
//output: dataframe summary {caption: Cultivation summary}
export function summary(summaryInput: DG.DataFrame): DG.DataFrame {
  const result = summaryInput.clone();
  result.name = 'Cultivation summary';
  return result;
}

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
