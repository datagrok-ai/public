import {getOptimizer, listOptimizers} from '..';
import {sphere} from './helpers';

describe('Registry', () => {
  it('lists registered optimizers', () => {
    const names = listOptimizers();
    expect(names).toContain('nelder-mead');
    expect(names).toContain('pso');
  });

  it('retrieves optimizer by name (case-insensitive)', () => {
    const opt = getOptimizer('Nelder-Mead');
    const r = opt.minimize(sphere, new Float64Array([5, -3, 7]), {
      maxIterations: 5_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(0, 4);
  });

  it('throws on unknown optimizer', () => {
    expect(() => getOptimizer('unknown')).toThrow('Unknown optimizer');
  });
});
