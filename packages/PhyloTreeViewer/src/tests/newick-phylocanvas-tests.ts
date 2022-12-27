import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';
import {parseNewick} from '@datagrok-libraries/bio';

category('newickParser_phylocanvas', () => {
  const nwk0 = `;`;
  const nwk1NoNameNoHeight = '();';
  const nwk1NameNoHeight = '(single);';
  const nwk1NameHeight = '(single:1.2);';
  const nwk1NoNameHeight = '(:1.2);';

  const nwk3LeafsNoHeight = `((n1,n2),n3);`;
  const nwk3LeafsIntNodesNoHeight = `((n1,n2)in1,n3);`;

  const enum Tests {
    nwk0 = 'nwk0',
    nwk1NoNameNoHeight = 'nwk1NoNameNoHeight',
    nwk1NameNoHeight = 'nwk1NameNoHeight',
    nwk1NameHeight = 'nwk1NameHeight',
    nwk1NoNameHeight = 'nwk1NoNameHeight',
    nwk3LeafsNoHeight = 'nwk3LeafsNoHeight',
    nwk3LeafsIntNodesNoHeight = 'nwk3LeafsIntNodesNoHeight',
  }

  const data: { [test: string]: { nwk: string, obj: Object } } = {
    [Tests.nwk0]: {
      nwk: ';',
      obj: {},
    },
    [Tests.nwk1NoNameNoHeight]: {
      nwk: '();',
      obj: {
        name: '', // root
        children: [{name: ''}]
      },
    },
    [Tests.nwk1NameNoHeight]: {
      nwk: '(single);',
      obj: {
        name: '', // root
        children: [{name: 'single'}]
      },
    },
    [Tests.nwk1NameHeight]: {
      nwk: '(single:1.2);',
      obj: {
        name: '', // root
        children: [{name: 'single', branch_length: 1.2}]
      },
    },
    [Tests.nwk1NoNameHeight]: {
      nwk: '(:1.2);',
      obj: {
        name: '', // root
        children: [{name: '', branch_length: 1.2}]
      },
    },
    [Tests.nwk3LeafsNoHeight]: {
      nwk: '((n1,n2),n3);',
      obj: {
        name: '', // root
        children: [
          {
            name: '',
            children: [{name: 'n1'}, {name: 'n2'}],
          },
          {name: 'n3'},
        ]
      },
    },
    [Tests.nwk3LeafsIntNodesNoHeight]: {
      nwk: '((n1,n2)in1,n3);',
      obj: {
        name: '', // root
        children: [
          {
            name: 'in1',
            children: [{name: 'n1'}, {name: 'n2'}],
          },
          {name: 'n3'},
        ]
      },
    },
  };

  test('nwk0', async () => {
    const testData = data[Tests.nwk0];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk1NoNameNoHeight', async () => {
    // const res = Newick.parse_newick(nwk1NoNameNoHeight);
    // let k = 11;
    // expectObject(res, {
    //   name: '', // root
    //   children: [{name: '',}]
    // });
    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk1NameNoHeight', async () => {
    // const res = Newick.parse_newick(nwk1NameNoHeight);
    // expectObject(res, {
    //   name: '', // root
    //   children: [{name: 'single',}]
    // });

    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk1NameHeight', async () => {
    // const res = Newick.parse_newick(nwk1NameHeight);
    // expectObject(res, {
    //   name: '', // root
    //   children: [{name: 'single', branch_length: 1.2}]
    // });

    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk1NoNameHeight', async () => {
    // const res = Newick.parse_newick(nwk1NoNameHeight);
    // expectObject(res, {
    //   name: '', // root
    //   children: [{name: '', branch_length: 1.2},]
    // });
    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk3LeafsNoHeight', async () => {
    // const res = Newick.parse_newick(nwk3LeafsNoHeight);
    // expectObject(res, {
    //   name: '', // root
    //   children: [
    //     {
    //       name: '',
    //       children: [{name: 'n1'}, {name: 'n2'},],
    //     },
    //     {name: 'n3',},
    //   ]
    // });

    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  test('nwk3LeafsIntNodesNoHeight', async () => {
    // const res = Newick.parse_newick(nwk3LeafsIntNodesNoHeight);
    // expectObject(res, {
    //   name: '', // root
    //   children: [
    //     {
    //       name: 'in1',
    //       children: [{name: 'n1'}, {name: 'n2'},],
    //     },
    //     {name: 'n3',},
    //   ]
    // });
    const testData = data[Tests.nwk1NoNameNoHeight];
    _testNewickToObject(testData.nwk, testData.obj);
  });

  function _testNewickToObject(nwk: string, tgtObj: Object) {
    const resObj = parseNewick(nwk);
    expectObject(resObj, tgtObj);
  }

  // function _testNodeToNewick(node:Object, tgtNwk: string){
  //   const resNwk = Newick.parse_newick()
  // }
});
