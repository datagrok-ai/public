import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {simulateDay} from '../bioreactor-model';

const closeTo = (actual: number, expected: number, tolerance = 1e-6): void => {
  expect(
    Math.abs(actual - expected) <= tolerance,
    true,
    `Expected ${actual} to be within ${tolerance} of ${expected}`,
  );
};

category('Bioreactor day calculation', () => {
  test('applies feeding and dilution before cultivation', async () => {
    const result = simulateDay(1, 3, 1, 2, 4, 0.4, 0.1, 0.5, 0.01, 0);

    closeTo(result.finalVolume, 2);
    closeTo(result.finalBiomass, 1);
    closeTo(result.finalSubstrate, 3.5);
    closeTo(result.finalBiomassMass, 2);
    closeTo(result.finalSubstrateMass, 7);
    expect(result.profile.length, 1);
  });

  test('decays biomass after substrate depletion', async () => {
    const result = simulateDay(0, 0, 1, 2, 0, 0.4, 0.1, 0.5, 0.02, 24);

    closeTo(result.finalBiomass, 2 * Math.exp(-0.02 * 24));
    closeTo(result.finalSubstrate, 0);
    closeTo(result.substrateConsumed, 0);
    expect(result.profile.length, 25);
  });

  test('preserves the growth yield mass balance without decay', async () => {
    const yieldCoefficient = 0.5;
    const result = simulateDay(0.25, 2, 1, 0.1, 0.5, 0.4, 0.1, yieldCoefficient, 0, 24);

    closeTo(result.biomassMassChange, yieldCoefficient * result.substrateConsumed, 1e-5);
    expect(result.finalBiomass >= 0, true);
    expect(result.finalSubstrate >= 0, true);
    expect(result.profile.some((point) => point.time % 1 !== 0 && point.substrate === 0), true);
  });
});
