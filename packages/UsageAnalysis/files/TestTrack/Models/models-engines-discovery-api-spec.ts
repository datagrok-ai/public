import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Models — Engine-discovery JS API (initEngines / getEngine / getAll)', async ({page}) => {
  test.setTimeout(180_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await softStep('S1: discovery scan returns non-empty + element-wise idempotent', async () => {
    const r = await page.evaluate(async () => {
      const DG = (window as any).DG;
      const scan = () => DG.Func.find({})
        .filter((f: any) => f.options && f.options.mlrole === 'train')
        .map((f: any) => f.options.mlname)
        .sort();
      const first = scan();
      const second = scan();
      return {
        firstCount: first.length,
        secondCount: second.length,
        firstNames: first,
        secondNames: second,
        equalElementwise: JSON.stringify(first) === JSON.stringify(second),
      };
    });
    expect(r.firstCount).toBeGreaterThan(0);
    expect(r.secondCount).toBe(r.firstCount);
    expect(r.equalElementwise).toBe(true);
    expect(r.firstNames).toEqual(r.secondNames);
  });

  await softStep('S2: at least one package-tagged train fn + matching apply sibling by mlname', async () => {
    const r = await page.evaluate(async () => {
      const DG = (window as any).DG;
      const all = DG.Func.find({});
      const train = all.filter((f: any) => f.options && f.options.mlrole === 'train');
      const apply = all.filter((f: any) => f.options && f.options.mlrole === 'apply');

      const trainByName: Record<string, any> = {};
      for (const f of train) trainByName[f.options.mlname] = {
        name: f.name, package: f.package && f.package.name, mlname: f.options.mlname,
      };
      const applyByName: Record<string, any> = {};
      for (const f of apply) applyByName[f.options.mlname] = {
        name: f.name, package: f.package && f.package.name, mlname: f.options.mlname,
      };

      const trainMlnames = Object.keys(trainByName);
      const pairedMlnames = trainMlnames.filter((n) => n in applyByName);

      const trainPkgs = new Set(train
        .map((f: any) => f.package && f.package.name)
        .filter((p: string | null | undefined) => !!p));

      return {
        trainCount: train.length,
        applyCount: apply.length,
        trainMlnames,
        pairedCount: pairedMlnames.length,
        unpairedTrainMlnames: trainMlnames.filter((n) => !(n in applyByName)),
        trainPackagesObserved: Array.from(trainPkgs).sort(),
      };
    });

    
    expect(r.trainCount).toBeGreaterThan(0);
    expect(r.pairedCount).toBeGreaterThan(0);

    
    expect(r.trainPackagesObserved.length).toBeGreaterThan(0);
  });

  
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
