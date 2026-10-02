import {category, test, expectArray} from '@datagrok-libraries/test/src/test';
import {parsePml, runPml} from '../pml-script';
import {CommandMessage, PatinaeViewer} from '../patinae-types';

function fakeViewer(calls: string[]): PatinaeViewer {
  return {
    init: async () => {},
    loadData: (data, name, format) => {calls.push(`load:${name}:${format}:${data.length}`);},
    execute: (command) => {calls.push(`exec:${command}`); return {messages: [{level: 'info', text: command}]};},
    executeAsync: async (command) => {
      calls.push(`exec:${command}`); return {messages: [{level: 'info', text: command}]};
    },
    getObjectInfos: () => [],
    destroy: () => {},
  };
}

category('Patinae: pml script', () => {
  test('parse', async () => {
    const lines = parsePml([
      '# comment', '', 'load 1crn.pdb, crambin', 'show cartoon, \\', '  crambin', '@more.pml',
      'load 4hhb.cif', 'load https://files.rcsb.org/download/1CRN.pdb', 'color red, all  ',
    ].join('\r\n'));
    expectArray(lines, [
      {kind: 'load', path: '1crn.pdb', name: 'crambin'},
      {kind: 'command', text: 'show cartoon, crambin'},
      {kind: 'include', path: 'more.pml'},
      {kind: 'load', path: '4hhb.cif', name: undefined},
      {kind: 'command', text: 'load https://files.rcsb.org/download/1CRN.pdb'},
      {kind: 'command', text: 'color red, all'},
    ]);
  });

  test('run resolves paths and includes', async () => {
    const calls: string[] = [];
    const files = {
      readText: async (path: string) => {calls.push(`text:${path}`); return 'hide everything';},
      readBytes: async (path: string) => {calls.push(`bytes:${path}`); return new Uint8Array(3);},
    };
    const messages: CommandMessage[] =
      await runPml('load a.pdb\n@inc.pml\nload b.cif.gz, bee\nshow sticks', 'System:AppData/Patinae/samples',
        files, fakeViewer(calls));
    expectArray(calls, [
      'bytes:System:AppData/Patinae/samples/a.pdb', 'load:a:pdb:3',
      'text:System:AppData/Patinae/samples/inc.pml', 'exec:hide everything',
      'bytes:System:AppData/Patinae/samples/b.cif.gz', 'load:bee:cif:3',
      'exec:show sticks',
    ]);
    expectArray(messages.map((m) => m.text), ['hide everything', 'show sticks']);
  });
});
