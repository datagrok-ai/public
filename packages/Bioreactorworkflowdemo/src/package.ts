/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import type {PipelineConfiguration} from '@datagrok-libraries/compute-api';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';
import timezone from 'dayjs/plugin/timezone';
export * from './package.g';

dayjs.extend(utc);
dayjs.extend(timezone);

export const _package = new DG.Package();

//name: BioreactorWorkflow
//description: Configure and simulate daily bacterial growth in a bioreactor
//meta.role: model
//editor: Compute2:TreeWizardEditor
//input: object params
//output: object result
export function bioreactorWorkflow(params: {version?: string}): PipelineConfiguration {
  void params;

  return {
    id: 'bioreactor-workflow',
    nqName: 'BioreactorWorkflowDemo:BioreactorWorkflow',
    version: '1.0',
    type: 'static',
    friendlyName: 'Bioreactor workflow demo',
    steps: [{
      id: 'bioreactor-configuration',
      nqName: 'BioreactorWorkflowDemo:BioreactorConfiguration',
      friendlyName: 'Bioreactor configuration',
    }],
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

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
