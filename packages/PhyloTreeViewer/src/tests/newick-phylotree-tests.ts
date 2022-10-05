import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import newickParser from 'phylotree/src/formats/newick';

category('newickParser_phylotree', () => {

  const nwk0 = `;`;
  const nwk1NoNameNoHeight = '();';
  const nwk1NameNoHeight = '(single);';
  const nwk1NameHeight = '(single:1.2);';
  const nwk1NoNameHeight = '(:1.2);';

  const nwk3LeafsNoHeight = `((n1,n2),n3);`;
  const nwk3LeafsIntNodesNoHeight = `((n1,n2)in1,n3);`;

  test('nwk0', async () => {
    const res = newickParser(nwk0);
    let k = 11;

    expect(res.error, null);
    expect(res.json.name, 'root');
    expect(res.json.children, undefined);
  });

  test('nwk1NoNameNoHeight', async () => {
    const res = newickParser(nwk1NoNameNoHeight);
    let k = 11;
    expectObject(res, {
      json: {
        name: 'root', children: [
          {name: '',}]
      },
      error: null,
    });
  });

  test('nwk1NameNoHeight', async () => {
    const res = newickParser(nwk1NameNoHeight);
    expectObject(res, {
      json: {
        name: 'root', children: [
          {name: 'single',},
        ]
      },
      error: null,
    });
  });

  test('nwk1NameHeight', async () => {
    const res = newickParser(nwk1NameHeight);
    expectObject(res, {
      json: {
        name: 'root', children: [
          {name: 'single', attribute: 1.2},
        ]
      },
      error: null,
    });
  });

  test('nwk1NoNameHeight', async () => {
    const res = newickParser(nwk1NoNameHeight);
    expectObject(res, {
      json: {
        name: 'root', children: [
          {name: '', attribute: 1.2},
        ]
      },
      error: null,
    });
  });

  test('nwk3LeafsNoHeight', async () => {
    const res = newickParser(nwk3LeafsNoHeight);
    expectObject(res, {
      json: {
        name: 'root', children: [
          {
            name: '',
            children: [{name: 'n1'}, {name: 'n2'},],
          },
          {name: 'n3',},
        ]
      },
      error: null,
    });
  });

  test('nwk3LeafsIntNodesNoHeight', async () => {
    const res = newickParser(nwk3LeafsIntNodesNoHeight);
    expectObject(res, {
      json: {
        name: 'root', children: [
          {
            name: 'in1',
            children: [{name: 'n1'}, {name: 'n2'},],
          },
          {name: 'n3',},
        ]
      },
      error: null,
    });
  });
});