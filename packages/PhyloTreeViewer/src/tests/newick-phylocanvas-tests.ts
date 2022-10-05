import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {Newick, NodeType} from '@phylocanvas/phylocanvas.gl';

category('newickParser_phylocanvas', () => {
  const nwk0 = `;`;
  const nwk1NoNameNoHeight = '();';
  const nwk1NameNoHeight = '(single);';
  const nwk1NameHeight = '(single:1.2);';
  const nwk1NoNameHeight = '(:1.2);';

  const nwk3LeafsNoHeight = `((n1,n2),n3);`;
  const nwk3LeafsIntNodesNoHeight = `((n1,n2)in1,n3);`;

  test('nwk0', async () => {
    const res = Newick.parse_newick(nwk0);
    expectObject(res, {});
  });

  test('nwk1NoNameNoHeight', async () => {
    const res = Newick.parse_newick(nwk1NoNameNoHeight);
    let k = 11;
    expectObject(res, {
      name: '', // root
      children: [{name: '',}]
    });
  });

  test('nwk1NameNoHeight', async () => {
    const res = Newick.parse_newick(nwk1NameNoHeight);
    expectObject(res, {
      name: '', // root
      children: [{name: 'single',}]
    });
  });

  test('nwk1NameHeight', async () => {
    const res = Newick.parse_newick(nwk1NameHeight);
    expectObject(res, {
      name: '', // root
      children: [{name: 'single', branch_length: 1.2}]
    });
  });

  test('nwk1NoNameHeight', async () => {
    const res = Newick.parse_newick(nwk1NoNameHeight);
    expectObject(res, {
      name: '', // root
      children: [{name: '', branch_length: 1.2},]
    });
  });

  test('nwk3LeafsNoHeight', async () => {
    const res = Newick.parse_newick(nwk3LeafsNoHeight);
    expectObject(res, {
      name: '', // root
      children: [
        {
          name: '',
          children: [{name: 'n1'}, {name: 'n2'},],
        },
        {name: 'n3',},
      ]
    });
  });

  test('nwk3LeafsIntNodesNoHeight', async () => {
    const res = Newick.parse_newick(nwk3LeafsIntNodesNoHeight);
    expectObject(res, {
      name: '', // root
      children: [
        {
          name: 'in1',
          children: [{name: 'n1'}, {name: 'n2'},],
        },
        {name: 'n3',},
      ]
    });
  });

});