import type {DojoWindowType} from './types';

declare const window: Window & DojoWindowType;

export async function initDojo(): Promise<void> {
  window.dojo$ = window.dojo$ || {};
  if (!window.dojo$.initPromise) {
    window.dojo$.initPromise = (async () => {
      const logPrefix: string = `dojo: _package.initDojo()`;
      console.debug(`${logPrefix}, start`);

      console.debug(`${logPrefix}, loadModules(), before`);
      const t1 = window.performance.now();
      await loadModules();
      const t2 = window.performance.now();
      console.debug(`${logPrefix}, loadModules(), after, ET: ${t2 - t1} ms`);

      console.debug(`${logPrefix}, end`);
    })();
  }
  return window.dojo$.initPromise;
}

const dojoPath = '../../vendor/dojo-1.10.10';

async function loadModules(): Promise<void> {
  const logPrefix: string = `dojo: _package.loadModules()`;

  const unc = '';
  window.dojo$.ctx = require.context('../../vendor/dojo-1.10.10', true, /\.js$/);
  // const unc = '.uncompressed.js'; // '';
  // window.dojo$.ctx = require.context('../../vendor/dojo-1.10.10', true, /\.js\.uncompressed\.js$/);

  window.dojo$.uncompressed = unc;
  await window.dojo$.ctx(`./dojo/dojo.js${unc}`);
}


(async () => {
  await initDojo();
})();
