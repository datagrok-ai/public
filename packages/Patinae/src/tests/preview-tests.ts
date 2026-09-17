import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {openPymol} from '../patinae-view';
import {_package} from '../package-test';

category('Patinae: preview', () => {
  for (const fn of ['1crn.prs', '1crn.pse', 'demo.pml', 'hiv-protease-tour.pml', 'trp-cage-ensemble.pml',
    'crambin-density.pml', 'protease-superposition.pml', 'hiv-protease.pse', 'trp-cage-ensemble.pse']) {
    test(`open ${fn}`, async () => {
      const file = (await grok.dapi.files.list(`System:AppData/${_package.name}/samples`, false, fn))[0];
      const {view, loaded} = openPymol(file);
      grok.shell.addView(view);
      try {
        const {viewer, messages} = await loaded;
        expect(viewer.getObjectInfos().length > 0, true, `${fn}: no objects loaded`);
        expect(messages.filter((m) => m.level === 'error').length, 0, `${fn}: script errors`);
      } finally {
        view.close();
      }
    }, 'gpu' in navigator ? {} : {skipReason: 'requires WebGPU'});
  }
});
