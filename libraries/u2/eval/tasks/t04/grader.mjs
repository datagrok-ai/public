import {controlOf, expect, flush, mount, need, offered} from '../../lib/harness.mjs';
import {frame} from '../../lib/platform.mjs';

const mol = (name) => ({name, type: 'string', semType: 'Molecule'});

const TABLES = {
  fresh: () => frame([{name: 'ID', type: 'string'}, mol('Molecule'), {name: 'MW', type: 'double'},
    mol('~Molecule.canonicalSmiles')],
  [{'ID': 'V1', 'Molecule': 'CCO', 'MW': 46, '~Molecule.canonicalSmiles': 'CCO'}], 'hits'),
  reloaded: () => frame([mol('~Core.canonicalSmiles'), {name: 'ID', type: 'string'}, mol('Core'), mol('Scaffold')],
    [{'~Core.canonicalSmiles': 'c1ccccc1', 'ID': 'V2', 'Core': 'c1ccccc1', 'Scaffold': 'C1CC1'}], 'campaign'),
};

async function picker(mod, table) {
  const input = need(mod, 'moleculeColumnInput')(TABLES[table]());
  const root = mount(input);
  await flush();
  const control = controlOf(root);
  expect(control !== undefined, 'the input is not a u2 control');
  return {value: control.getState().value, items: offered(control).map((i) => i.key)};
}

export const checks = {
  'offers-molecules': async (mod) => {
    for (const [table, wanted] of [['fresh', ['Molecule']], ['reloaded', ['Core', 'Scaffold']]]) {
      const {items} = await picker(mod, table);
      for (const name of wanted)
        expect(items.includes(name), `the ${table} table's molecule column ${name} is not among ${items.join(', ')}`);
      expect(!items.includes('ID'), 'a column that is not a molecule is offered');
    }
  },
  'no-hidden-columns': async (mod) => {
    for (const table of Object.keys(TABLES)) {
      const {items} = await picker(mod, table);
      const hidden = items.filter((n) => String(n).startsWith('~'));
      expect(hidden.length === 0, `offers ${hidden.join(', ')}, a column the platform keeps out of every list`);
    }
  },
  'preselects-molecule': async (mod) => {
    for (const [table, wanted] of [['fresh', 'Molecule'], ['reloaded', 'Core']]) {
      const {value} = await picker(mod, table);
      expect(value === wanted, `on the ${table} table it starts on ${JSON.stringify(value)}, not ${wanted}`);
    }
  },
};
