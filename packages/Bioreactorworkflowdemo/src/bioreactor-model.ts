export type BioreactorState = {
  biomass: number;
  substrate: number;
};

export type DailyProfilePoint = {
  time: number;
  volume: number;
  biomass: number;
  substrate: number;
  biomassMass: number;
  substrateMass: number;
};

export type DaySimulationResult = {
  finalVolume: number;
  finalBiomass: number;
  finalSubstrate: number;
  finalBiomassMass: number;
  finalSubstrateMass: number;
  biomassMassChange: number;
  substrateConsumed: number;
  profile: DailyProfilePoint[];
};

const SOLVER_TOLERANCE = 1e-8;
const MIN_STEP = 1e-8;
const MAX_STEP = 0.25;

type IntegrationResult = {
  state: BioreactorState;
  depletion?: {
    time: number;
    state: BioreactorState;
  };
};

function derivatives(
  state: BioreactorState,
  maximumGrowthRate: number,
  halfSaturationConstant: number,
  yieldCoefficient: number,
  decayConstant: number,
): BioreactorState {
  const biomass = Math.max(0, state.biomass);
  const substrate = Math.max(0, state.substrate);
  const growthRate = substrate === 0 ? 0 :
    maximumGrowthRate * substrate / (halfSaturationConstant + substrate);

  return {
    biomass: (growthRate - decayConstant) * biomass,
    substrate: -growthRate * biomass / yieldCoefficient,
  };
}

function rk4Step(
  state: BioreactorState,
  step: number,
  maximumGrowthRate: number,
  halfSaturationConstant: number,
  yieldCoefficient: number,
  decayConstant: number,
): BioreactorState {
  const derivative = (value: BioreactorState) => derivatives(
    value,
    maximumGrowthRate,
    halfSaturationConstant,
    yieldCoefficient,
    decayConstant,
  );
  const addScaled = (value: BioreactorState, slope: BioreactorState, scale: number): BioreactorState => ({
    biomass: value.biomass + slope.biomass * scale,
    substrate: value.substrate + slope.substrate * scale,
  });

  const k1 = derivative(state);
  const k2 = derivative(addScaled(state, k1, step / 2));
  const k3 = derivative(addScaled(state, k2, step / 2));
  const k4 = derivative(addScaled(state, k3, step));

  return {
    biomass: state.biomass + step * (k1.biomass + 2 * k2.biomass + 2 * k3.biomass + k4.biomass) / 6,
    substrate: state.substrate + step *
      (k1.substrate + 2 * k2.substrate + 2 * k3.substrate + k4.substrate) / 6,
  };
}

function integrateTo(
  initialState: BioreactorState,
  startTime: number,
  endTime: number,
  maximumGrowthRate: number,
  halfSaturationConstant: number,
  yieldCoefficient: number,
  decayConstant: number,
): IntegrationResult {
  let state = {...initialState};
  let time = startTime;
  let step = Math.min(MAX_STEP, endTime - startTime);
  let depletion: IntegrationResult['depletion'];

  while (time < endTime) {
    step = Math.min(step, endTime - time);

    if (state.substrate <= SOLVER_TOLERANCE) {
      state = {
        biomass: state.biomass * Math.exp(-decayConstant * (endTime - time)),
        substrate: 0,
      };
      break;
    }

    const fullStep = rk4Step(
      state, step, maximumGrowthRate, halfSaturationConstant, yieldCoefficient, decayConstant,
    );
    const firstHalf = rk4Step(
      state, step / 2, maximumGrowthRate, halfSaturationConstant, yieldCoefficient, decayConstant,
    );
    const twoHalfSteps = rk4Step(
      firstHalf, step / 2, maximumGrowthRate, halfSaturationConstant, yieldCoefficient, decayConstant,
    );
    const biomassScale = SOLVER_TOLERANCE * Math.max(1, Math.abs(state.biomass), Math.abs(twoHalfSteps.biomass));
    const substrateScale = SOLVER_TOLERANCE *
      Math.max(1, Math.abs(state.substrate), Math.abs(twoHalfSteps.substrate));
    const error = Math.max(
      Math.abs(twoHalfSteps.biomass - fullStep.biomass) / (15 * biomassScale),
      Math.abs(twoHalfSteps.substrate - fullStep.substrate) / (15 * substrateScale),
    );

    if (error <= 1 || step <= MIN_STEP) {
      if (twoHalfSteps.substrate <= SOLVER_TOLERANCE) {
        let lower = 0;
        let upper = step;
        for (let i = 0; i < 40; i++) {
          const middle = (lower + upper) / 2;
          const middleState = rk4Step(
            state, middle, maximumGrowthRate, halfSaturationConstant, yieldCoefficient, decayConstant,
          );
          if (middleState.substrate > SOLVER_TOLERANCE)
            lower = middle;
          else
            upper = middle;
        }
        const depletionState = rk4Step(
          state, upper, maximumGrowthRate, halfSaturationConstant, yieldCoefficient, decayConstant,
        );
        depletion = {
          time: time + upper,
          state: {biomass: Math.max(0, depletionState.biomass), substrate: 0},
        };
        state = {
          biomass: depletion.state.biomass * Math.exp(-decayConstant * (endTime - depletion.time)),
          substrate: 0,
        };
        break;
      } else {
        state = {
          biomass: Math.max(0, twoHalfSteps.biomass),
          substrate: twoHalfSteps.substrate,
        };
        time += step;
      }
    }

    const factor = error === 0 ? 2 : Math.min(2, Math.max(0.2, 0.9 * Math.pow(error, -0.2)));
    step = Math.max(MIN_STEP, Math.min(MAX_STEP, step * factor));
  }

  return {state, depletion};
}

function validateInputs(values: Record<string, number>, positive: string[]): void {
  for (const [name, value] of Object.entries(values)) {
    if (!Number.isFinite(value))
      throw new Error(`${name} must be finite`);
    if (value < 0 || (positive.includes(name) && value === 0))
      throw new Error(`${name} must be ${positive.includes(name) ? 'positive' : 'non-negative'}`);
  }
}

export function simulateDay(
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
): DaySimulationResult {
  const values = {
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
  };
  validateInputs(values, [
    'incomingVolume',
    'maximumGrowthRate',
    'halfSaturationConstant',
    'yieldCoefficient',
  ]);

  const finalVolume = incomingVolume + solutionAdded;
  const incomingBiomassMass = incomingBiomass * incomingVolume;
  const availableSubstrateMass = incomingSubstrate * incomingVolume + substrateAdded;
  let state: BioreactorState = {
    biomass: incomingBiomassMass / finalVolume,
    substrate: availableSubstrateMass / finalVolume,
  };
  const profile: DailyProfilePoint[] = [];
  const addProfilePoint = (time: number) => profile.push({
    time,
    volume: finalVolume,
    biomass: state.biomass,
    substrate: state.substrate,
    biomassMass: state.biomass * finalVolume,
    substrateMass: state.substrate * finalVolume,
  });

  addProfilePoint(0);
  let time = 0;
  while (time < dayDuration) {
    const nextTime = Math.min(dayDuration, Math.floor(time) + 1);
    const integration = integrateTo(
      state,
      time,
      nextTime,
      maximumGrowthRate,
      halfSaturationConstant,
      yieldCoefficient,
      decayConstant,
    );
    if (integration.depletion != null) {
      state = integration.depletion.state;
      addProfilePoint(integration.depletion.time);
    }
    state = integration.state;
    time = nextTime;
    addProfilePoint(time);
  }

  const finalBiomassMass = state.biomass * finalVolume;
  const finalSubstrateMass = state.substrate * finalVolume;
  return {
    finalVolume,
    finalBiomass: state.biomass,
    finalSubstrate: state.substrate,
    finalBiomassMass,
    finalSubstrateMass,
    biomassMassChange: finalBiomassMass - incomingBiomassMass,
    substrateConsumed: availableSubstrateMass - finalSubstrateMass,
    profile,
  };
}
