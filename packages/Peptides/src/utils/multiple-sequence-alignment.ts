/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/bio/src/sequence-encoder';
import {FSInjector} from './fs-injector';

const initKalign = require('../wasm/kalign');

export async function runKalign(col: DG.Column, webRootValue: string) : Promise<DG.Column> {
  const _webPath = `${webRootValue}dist/kalign.wasm`;
  const sequences = col.toList().map((v: string, _) => AlignedSequenceEncoder.clean(v).replace(/\-/g, ''));
  const stdin = sequences.reduce((a, v, i) => a + `>sample${i}\n${v}\n`);
  const inj = new FSInjector(stdin);

  console.log(initKalign);
  console.log(_webPath);

  let kalign = await initKalign({locateFile: () => _webPath});

  const Module = {
    preRun: function() {
      console.log('preRun');
      inj.reset();
      kalign.FS.init(inj.stdin.bind(inj), inj.stdout.bind(inj), inj.stderr.bind(inj));
    },
    locateFile: () => _webPath,
  };

  console.log(kalign);
  console.log(Module);

  kalign = await initKalign(Module);
  console.log(kalign);
  await kalign.run();
  console.log(inj);

  return DG.Column.fromStrings('aligned', inj.stdoutBuffer.split(/>sample\d+\n/g));
}
