export const SAMPLE_TREE = [
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[c:1]([cH2:2])([cH2:3])[C:5]([CH3:4])[CH:6]=[O:7]>>Cl[c:1]([cH2:2])[cH2:3].[CH3:4][C:5][CH:6]=[O:7]',
        'is_reaction': true,
        'metadata': {
          'template_hash': 'fb6366c2f1f5cc39b7a6b203cef6ff4bbb326cc072fa4ca9c0153f198820a553',
          'classification': '0.0 Unrecognized',
          'library_occurence': 30,
          'policy_probability': 0.0009,
          'policy_probability_rank': 17,
          'policy_name': 'uspto',
          'template_code': 41811,
          'template': '[C:4]-[CH;D3;+0:5](-[C:6]=[O;D1;H0:7])-[c;H0;D3;+0:1](:[c:2]):[c:3]>>Cl-[c;H0;D3;+0:1](:[c:2]):[c:3].[C:4]-[CH2;D2;+0:5]-[C:6]=[O;D1;H0:7]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH2:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH3:22].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([Cl:27])[cH:23][cH:24][c:25]2[cH:26]1',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(Cl)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[C:1][O:2][cH3:3]>>COS(=O)(=O)O[C:1].[O:2][cH3:3]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '7572d81cb8fe1802366cd5108e7ad9313d2ce90888f05fbdc714a1b2fc719786',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 511,
                  'policy_probability': 0.1630000025,
                  'policy_probability_rank': 0,
                  'policy_name': 'uspto',
                  'template_code': 19597,
                  'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-O-S(=O)(=O)-O-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([Cl:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:9][S:10]([O:11][CH3:12])(=[O:13])=[O:14].[OH:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([Cl:27])[cH:23][cH:24][c:25]2[cH:26]1',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COS(=O)(=O)OC',
                    'is_chemical': true,
                    'in_stock': true,
                  },
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'Oc1ccc2cc(Cl)ccc2c1',
                    'is_chemical': true,
                    'in_stock': true,
                  },
                ],
              },
            ],
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'CCC(=O)Oc1ccc(C)cc1OC',
            'is_chemical': true,
            'in_stock': true,
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.9940398539,
      'number of reactions': 2,
      'number of pre-cursors': 3,
      'number of pre-cursors in stock': 3,
      'average template occurrence': 270.5,
    },
    'metadata': {
      'created_at_iteration': 11,
      'is_solved': true,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>COS(=O)(=O)O[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '7572d81cb8fe1802366cd5108e7ad9313d2ce90888f05fbdc714a1b2fc719786',
          'classification': '0.0 Unrecognized',
          'library_occurence': 511,
          'policy_probability': 0.0046000001,
          'policy_probability_rank': 6,
          'policy_name': 'uspto',
          'template_code': 19597,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-O-S(=O)(=O)-O-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][O:27][S:28]([O:29][CH3:30])(=[O:31])=[O:32]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COS(=O)(=O)OC',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[cH3:1][O:2]>>c1ccc(C[O:2][cH3:1])cc1',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '69a58578a3918c6cc67c5d5ac4a9bb4396706cbfeeccebc3be07c7134dbfc78f',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 6005,
                  'policy_probability': 0.0322999991,
                  'policy_probability_rank': 3,
                  'policy_name': 'uspto',
                  'template_code': 17606,
                  'template': '[OH;D1;+0:2]-[c:1]>>[c:1]-[O;H0;D2;+0:2]-C-c1:c:c:c:c:c:1',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH2:27][c:28]3[cH:29][cH:30][cH:31][cH:32][cH:33]3)[cH:23][cH:24][c:25]2[cH:26]1',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OCc3ccccc3)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[c:1]([cH2:2])([cH2:3])[O:6][CH:5]=[O:4]>>O[c:1]([cH2:2])[cH2:3].[O:4]=[CH:5][O:6]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '4af293754194d4293822dac15dffac886aafe5cce7291212bae7605614f9968a',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 737,
                          'policy_probability': 0.2101999968,
                          'policy_probability_rank': 1,
                          'policy_name': 'uspto',
                          'template_code': 12567,
                          'template': '[O;D1;H0:4]=[C:5]-[O;H0;D2;+0:6]-[c;H0;D3;+0:1](:[c:2]):[c:3]>>O-[c;H0;D3;+0:1](:[c:2]):[c:3].[O;D1;H0:4]=[C:5]-[OH;D1;+0:6]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH2:27][c:28]3[cH:29][cH:30][cH:31][cH:32][cH:33]3)[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:13])[cH:23][cH:24][c:25]2[cH:26]1.[c:14]1([OH:34])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'Cc1ccc(O)c(OCc2ccccc2)c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][cH3:2]>>CC(=O)[O:1][cH3:2]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': 'ba9cf922694c22786e94d6562607c5b880cd83d9f9ace2090b7b33f921f84f06',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 964,
                                  'policy_probability': 0.2029000074,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 31122,
                                  'template': '[OH;D1;+0:1]-[c:2]>>C-C(=O)-[O;H0;D2;+0:1]-[c:2]',
                                  'mapped_reaction_smiles': '[c:14]1([OH:34])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[c:14]1([O:34][C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'CC(=O)Oc1ccc(C)cc1OCc1ccccc1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[O:1]([C:3]([CH3:2])=[O:4])[c:5]([cH2:6])[cH2:7]>>O=C(O[O:1])c1cccc(Cl)c1.[CH3:2][C:3](=[O:4])[c:5]([cH2:6])[cH2:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '1ec0ea5b4cb26ec0598139be26b1ffece0af34c9422e1d0ac10d34892739d5c8',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 104,
                                          'policy_probability': 0.0644999966,
                                          'policy_probability_rank': 2,
                                          'policy_name': 'uspto',
                                          'template_code': 5183,
                                          'template': '[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[O;H0;D2;+0:1]-[c;H0;D3;+0:5](:[c:6]):[c:7]>>Cl-c1:c:c:c:c(-C(=O)-O-[OH;D1;+0:1]):c:1.[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[c;H0;D3;+0:5](:[c:6]):[c:7]',
                                          'mapped_reaction_smiles': '[c:14]1([O:34][C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[OH:34][O:38][C:39]([c:40]1[cH:41][cH:42][cH:43][c:44]([Cl:45])[cH:46]1)=[O:47].[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'O=C(OO)c1cccc(Cl)c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CC(=O)c1ccc(C)cc1OCc1ccccc1',
                                            'is_chemical': true,
                                            'in_stock': false,
                                            'children': [
                                              {
                                                'type': 'reaction',
                                                'hide': false,
                                                'smiles': '[C:1]([cH3:2])[O:3][cH3:4]>>Br[C:1][cH3:2].[O:3][cH3:4]',
                                                'is_reaction': true,
                                                'metadata': {
                                                  'template_hash': '569f8f49ca4d8c992398ac28a432fca74adaf7a78d6a701f24550348c4be4c40',
                                                  'classification': '0.0 Unrecognized',
                                                  'library_occurence': 5146,
                                                  'policy_probability': 0.5411000252,
                                                  'policy_probability_rank': 0,
                                                  'policy_name': 'uspto',
                                                  'template_code': 14455,
                                                  'template': '[c:4]-[O;H0;D2;+0:3]-[CH2;D2;+0:1]-[c:2]>>Br-[CH2;D2;+0:1]-[c:2].[OH;D1;+0:3]-[c:4]',
                                                  'mapped_reaction_smiles': '[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[CH2:27]([c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1)[Br:34].[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[OH:21]',
                                                },
                                                'children': [
                                                  {
                                                    'type': 'mol',
                                                    'hide': false,
                                                    'smiles': 'CC(=O)c1ccc(C)cc1O',
                                                    'is_chemical': true,
                                                    'in_stock': true,
                                                  },
                                                  {
                                                    'type': 'mol',
                                                    'hide': false,
                                                    'smiles': 'BrCc1ccccc1',
                                                    'is_chemical': true,
                                                    'in_stock': true,
                                                  },
                                                ],
                                              },
                                            ],
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][C:2]([CH3:3])=[O:4]>>C[O:1][C:2]([CH3:3])=[O:4]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '110b771d77121b19f0a012d7c67263575768c8985321b84bf1629ea79c879d51',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 21925,
                                  'policy_probability': 0.4129999876,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 2818,
                                  'template': '[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:13])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COC(=O)C(C)c1ccc2cc(OC)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[C:1][C:6]([C:4]([O:3][CH3:2])=[O:5])[cH3:7]>>I[C:1].[CH3:2][O:3][C:4](=[O:5])[C:6][cH3:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '3c075bd2dba9541386d3ade4d2fb72e3c039b52e5142cf25a42d6cba6563011d',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 148,
                                          'policy_probability': 0.5683000088,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 10055,
                                          'template': '[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH;D3;+0:6](-[CH3;D1;+0:1])-[c:7]>>I-[CH3;D1;+0:1].[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH2;D2;+0:6]-[c:7]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:10][I:14].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CI',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COC(=O)Cc1ccc2cc(OC)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7976268128,
      'number of reactions': 8,
      'number of pre-cursors': 6,
      'number of pre-cursors in stock': 5,
      'average template occurrence': 4442.5,
    },
    'metadata': {
      'created_at_iteration': 74,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>CC(=O)[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '9acc680c930290c191eb4a208e5cf64573c699d51b010ffd4ddde6d301764851',
          'classification': '0.0 Unrecognized',
          'library_occurence': 19,
          'policy_probability': 0.0004,
          'policy_probability_rank': 29,
          'policy_name': 'uspto',
          'template_code': 25786,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-C(=O)-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][C:27]([CH3:28])=[O:29]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'CC(C)=O',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[cH3:1][O:2]>>c1ccc(C[O:2][cH3:1])cc1',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '69a58578a3918c6cc67c5d5ac4a9bb4396706cbfeeccebc3be07c7134dbfc78f',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 6005,
                  'policy_probability': 0.0322999991,
                  'policy_probability_rank': 3,
                  'policy_name': 'uspto',
                  'template_code': 17606,
                  'template': '[OH;D1;+0:2]-[c:1]>>[c:1]-[O;H0;D2;+0:2]-C-c1:c:c:c:c:c:1',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH2:27][c:28]3[cH:29][cH:30][cH:31][cH:32][cH:33]3)[cH:23][cH:24][c:25]2[cH:26]1',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OCc3ccccc3)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[c:1]([cH2:2])([cH2:3])[O:6][CH:5]=[O:4]>>O[c:1]([cH2:2])[cH2:3].[O:4]=[CH:5][O:6]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '4af293754194d4293822dac15dffac886aafe5cce7291212bae7605614f9968a',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 737,
                          'policy_probability': 0.2101999968,
                          'policy_probability_rank': 1,
                          'policy_name': 'uspto',
                          'template_code': 12567,
                          'template': '[O;D1;H0:4]=[C:5]-[O;H0;D2;+0:6]-[c;H0;D3;+0:1](:[c:2]):[c:3]>>O-[c;H0;D3;+0:1](:[c:2]):[c:3].[O;D1;H0:4]=[C:5]-[OH;D1;+0:6]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH2:27][c:28]3[cH:29][cH:30][cH:31][cH:32][cH:33]3)[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:13])[cH:23][cH:24][c:25]2[cH:26]1.[c:14]1([OH:34])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'Cc1ccc(O)c(OCc2ccccc2)c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][cH3:2]>>CC(=O)[O:1][cH3:2]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': 'ba9cf922694c22786e94d6562607c5b880cd83d9f9ace2090b7b33f921f84f06',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 964,
                                  'policy_probability': 0.2029000074,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 31122,
                                  'template': '[OH;D1;+0:1]-[c:2]>>C-C(=O)-[O;H0;D2;+0:1]-[c:2]',
                                  'mapped_reaction_smiles': '[c:14]1([OH:34])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[c:14]1([O:34][C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'CC(=O)Oc1ccc(C)cc1OCc1ccccc1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[O:1]([C:3]([CH3:2])=[O:4])[c:5]([cH2:6])[cH2:7]>>O=C(O[O:1])c1cccc(Cl)c1.[CH3:2][C:3](=[O:4])[c:5]([cH2:6])[cH2:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '1ec0ea5b4cb26ec0598139be26b1ffece0af34c9422e1d0ac10d34892739d5c8',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 104,
                                          'policy_probability': 0.0644999966,
                                          'policy_probability_rank': 2,
                                          'policy_name': 'uspto',
                                          'template_code': 5183,
                                          'template': '[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[O;H0;D2;+0:1]-[c;H0;D3;+0:5](:[c:6]):[c:7]>>Cl-c1:c:c:c:c(-C(=O)-O-[OH;D1;+0:1]):c:1.[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[c;H0;D3;+0:5](:[c:6]):[c:7]',
                                          'mapped_reaction_smiles': '[c:14]1([O:34][C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[OH:34][O:38][C:39]([c:40]1[cH:41][cH:42][cH:43][c:44]([Cl:45])[cH:46]1)=[O:47].[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'O=C(OO)c1cccc(Cl)c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CC(=O)c1ccc(C)cc1OCc1ccccc1',
                                            'is_chemical': true,
                                            'in_stock': false,
                                            'children': [
                                              {
                                                'type': 'reaction',
                                                'hide': false,
                                                'smiles': '[C:1]([cH3:2])[O:3][cH3:4]>>Br[C:1][cH3:2].[O:3][cH3:4]',
                                                'is_reaction': true,
                                                'metadata': {
                                                  'template_hash': '569f8f49ca4d8c992398ac28a432fca74adaf7a78d6a701f24550348c4be4c40',
                                                  'classification': '0.0 Unrecognized',
                                                  'library_occurence': 5146,
                                                  'policy_probability': 0.5411000252,
                                                  'policy_probability_rank': 0,
                                                  'policy_name': 'uspto',
                                                  'template_code': 14455,
                                                  'template': '[c:4]-[O;H0;D2;+0:3]-[CH2;D2;+0:1]-[c:2]>>Br-[CH2;D2;+0:1]-[c:2].[OH;D1;+0:3]-[c:4]',
                                                  'mapped_reaction_smiles': '[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[O:21][CH2:27][c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1>>[CH2:27]([c:28]1[cH:29][cH:30][cH:31][cH:32][cH:33]1)[Br:34].[c:14]1([C:35]([CH3:36])=[O:37])[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[OH:21]',
                                                },
                                                'children': [
                                                  {
                                                    'type': 'mol',
                                                    'hide': false,
                                                    'smiles': 'CC(=O)c1ccc(C)cc1O',
                                                    'is_chemical': true,
                                                    'in_stock': true,
                                                  },
                                                  {
                                                    'type': 'mol',
                                                    'hide': false,
                                                    'smiles': 'BrCc1ccccc1',
                                                    'is_chemical': true,
                                                    'in_stock': true,
                                                  },
                                                ],
                                              },
                                            ],
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][C:2]([CH3:3])=[O:4]>>C[O:1][C:2]([CH3:3])=[O:4]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '110b771d77121b19f0a012d7c67263575768c8985321b84bf1629ea79c879d51',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 21925,
                                  'policy_probability': 0.4129999876,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 2818,
                                  'template': '[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:13])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COC(=O)C(C)c1ccc2cc(OC)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[C:1][C:6]([C:4]([O:3][CH3:2])=[O:5])[cH3:7]>>I[C:1].[CH3:2][O:3][C:4](=[O:5])[C:6][cH3:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '3c075bd2dba9541386d3ade4d2fb72e3c039b52e5142cf25a42d6cba6563011d',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 148,
                                          'policy_probability': 0.5683000088,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 10055,
                                          'template': '[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH;D3;+0:6](-[CH3;D1;+0:1])-[c:7]>>I-[CH3;D1;+0:1].[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH2;D2;+0:6]-[c:7]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:10][I:14].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][CH3:27])[cH:23][cH:24][c:25]2[cH:26]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CI',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COC(=O)Cc1ccc2cc(OC)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7976268128,
      'number of reactions': 8,
      'number of pre-cursors': 6,
      'number of pre-cursors in stock': 5,
      'average template occurrence': 4381,
    },
    'metadata': {
      'created_at_iteration': 75,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][C:4]([CH:3]=[O:2])[cH3:5]>>I[C:1].[O:2]=[CH:3][C:4][cH3:5]',
        'is_reaction': true,
        'metadata': {
          'template_hash': 'e82f6c5effc3d594fdfcdda786c450919c8169e3007c596602b74e136b93d3eb',
          'classification': '0.0 Unrecognized',
          'library_occurence': 117,
          'policy_probability': 0.0003,
          'policy_probability_rank': 34,
          'policy_name': 'uspto',
          'template_code': 38578,
          'template': '[CH3;D1;+0:1]-[CH;D3;+0:4](-[c:5])-[C:3]=[O;D1;H0:2]>>I-[CH3;D1;+0:1].[O;D1;H0:2]=[C:3]-[CH2;D2;+0:4]-[c:5]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:10][I:27].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'CI',
            'is_chemical': true,
            'in_stock': false,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(CC(=O)Oc3ccc(C)cc3OC)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[C:1][O:2][cH3:3]>>COS(=O)(=O)O[C:1].[O:2][cH3:3]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '7572d81cb8fe1802366cd5108e7ad9313d2ce90888f05fbdc714a1b2fc719786',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 511,
                  'policy_probability': 0.0020999999,
                  'policy_probability_rank': 6,
                  'policy_name': 'uspto',
                  'template_code': 19597,
                  'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-O-S(=O)(=O)-O-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][O:27][S:28]([O:29][CH3:30])(=[O:31])=[O:32]',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COS(=O)(=O)OC',
                    'is_chemical': true,
                    'in_stock': true,
                  },
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(CC(=O)Oc3ccc(C)cc3O)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 1107,
                          'policy_probability': 0.6161000133,
                          'policy_probability_rank': 0,
                          'policy_name': 'uspto',
                          'template_code': 248,
                          'template': '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1.[OH:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[OH:21]',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'Cc1ccc(O)c(O)c1',
                            'is_chemical': true,
                            'in_stock': true,
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(CC(=O)Cl)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[Cl:1][C:2]([CH3:3])=[O:4]>>O=C(Cl)C(=O)[Cl:1].O[C:2]([CH3:3])=[O:4]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '866e1a8ac889f536ad4f364036065c6e47f9d15c5f533dbb707d5e47a3434b1a',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 1623,
                                  'policy_probability': 0.5088000298,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 22424,
                                  'template': '[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-C(=O)-C(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[OH:28])[cH:23][cH:24][c:25]2[cH:26]1.[Cl:27][C:29]([C:30]([Cl:31])=[O:32])=[O:33]',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'O=C(Cl)C(=O)Cl',
                                    'is_chemical': true,
                                    'in_stock': true,
                                  },
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COc1ccc2cc(CC(=O)O)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': true,
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.785,
      'number of reactions': 4,
      'number of pre-cursors': 5,
      'number of pre-cursors in stock': 4,
      'average template occurrence': 839.5,
    },
    'metadata': {
      'created_at_iteration': 83,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>COS(=O)(=O)O[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '7572d81cb8fe1802366cd5108e7ad9313d2ce90888f05fbdc714a1b2fc719786',
          'classification': '0.0 Unrecognized',
          'library_occurence': 511,
          'policy_probability': 0.0046000001,
          'policy_probability_rank': 6,
          'policy_name': 'uspto',
          'template_code': 19597,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-O-S(=O)(=O)-O-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][O:27][S:28]([O:29][CH3:30])(=[O:31])=[O:32]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COS(=O)(=O)OC',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 1107,
                  'policy_probability': 0.622600019,
                  'policy_probability_rank': 0,
                  'policy_name': 'uspto',
                  'template_code': 248,
                  'template': '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1.[OH:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[OH:21]',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'Cc1ccc(O)c(O)c1',
                    'is_chemical': true,
                    'in_stock': true,
                  },
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Cl)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 1863,
                          'policy_probability': 0.4481999874,
                          'policy_probability_rank': 0,
                          'policy_name': 'uspto',
                          'template_code': 4355,
                          'template': '[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:28])[cH:23][cH:24][c:25]2[cH:26]1.[Cl:27][S:29]([Cl:30])=[O:31]',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'O=S(Cl)Cl',
                            'is_chemical': true,
                            'in_stock': true,
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][C:2]([CH3:3])=[O:4]>>C[O:1][C:2]([CH3:3])=[O:4]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '110b771d77121b19f0a012d7c67263575768c8985321b84bf1629ea79c879d51',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 21925,
                                  'policy_probability': 0.4129999876,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 2818,
                                  'template': '[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:28])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COC(=O)C(C)c1ccc2cc(OC)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[C:1][C:6]([C:4]([O:3][CH3:2])=[O:5])[cH3:7]>>I[C:1].[CH3:2][O:3][C:4](=[O:5])[C:6][cH3:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '3c075bd2dba9541386d3ade4d2fb72e3c039b52e5142cf25a42d6cba6563011d',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 148,
                                          'policy_probability': 0.5683000088,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 10055,
                                          'template': '[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH;D3;+0:6](-[CH3;D1;+0:1])-[c:7]>>I-[CH3;D1;+0:1].[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH2;D2;+0:6]-[c:7]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:10][I:13].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CI',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COC(=O)Cc1ccc2cc(OC)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7734470711,
      'number of reactions': 5,
      'number of pre-cursors': 5,
      'number of pre-cursors in stock': 4,
      'average template occurrence': 5110.8,
    },
    'metadata': {
      'created_at_iteration': 4,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>COS(=O)(=O)O[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '7572d81cb8fe1802366cd5108e7ad9313d2ce90888f05fbdc714a1b2fc719786',
          'classification': '0.0 Unrecognized',
          'library_occurence': 511,
          'policy_probability': 0.0046000001,
          'policy_probability_rank': 6,
          'policy_name': 'uspto',
          'template_code': 19597,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-O-S(=O)(=O)-O-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][O:27][S:28]([O:29][CH3:30])(=[O:31])=[O:32]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COS(=O)(=O)OC',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[O:1][cH3:2]>>CC(=O)[O:1][cH3:2]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': 'ba9cf922694c22786e94d6562607c5b880cd83d9f9ace2090b7b33f921f84f06',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 964,
                  'policy_probability': 0.0196000002,
                  'policy_probability_rank': 4,
                  'policy_name': 'uspto',
                  'template_code': 31122,
                  'template': '[OH;D1;+0:1]-[c:2]>>C-C(=O)-[O;H0;D2;+0:1]-[c:2]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC(C)=O)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[O:1]([C:3]([CH3:2])=[O:4])[c:5]([cH2:6])[cH2:7]>>O=C(O[O:1])c1cccc(Cl)c1.[CH3:2][C:3](=[O:4])[c:5]([cH2:6])[cH2:7]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '1ec0ea5b4cb26ec0598139be26b1ffece0af34c9422e1d0ac10d34892739d5c8',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 104,
                          'policy_probability': 0.0810000002,
                          'policy_probability_rank': 3,
                          'policy_name': 'uspto',
                          'template_code': 5183,
                          'template': '[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[O;H0;D2;+0:1]-[c;H0;D3;+0:5](:[c:6]):[c:7]>>Cl-c1:c:c:c:c(-C(=O)-O-[OH;D1;+0:1]):c:1.[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[c;H0;D3;+0:5](:[c:6]):[c:7]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1.[OH:21][O:22][C:30]([c:31]1[cH:32][cH:33][cH:34][c:35]([Cl:36])[cH:37]1)=[O:38]',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'O=C(OO)c1cccc(Cl)c1',
                            'is_chemical': true,
                            'in_stock': true,
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3C(C)=O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 1107,
                                  'policy_probability': 0.6410999894,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 248,
                                  'template': '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:30])[cH:23][cH:24][c:25]2[cH:26]1.[OH:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[C:27]([CH3:28])=[O:29]',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'CC(=O)c1cc(C)ccc1O',
                                    'is_chemical': true,
                                    'in_stock': true,
                                  },
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COc1ccc2cc(C(C)C(=O)Cl)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 1863,
                                          'policy_probability': 0.4481999874,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 4355,
                                          'template': '[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:30])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:27])[cH:23][cH:24][c:25]2[cH:26]1.[Cl:30][S:31]([Cl:32])=[O:33]',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'O=S(Cl)Cl',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7734470711,
      'number of reactions': 5,
      'number of pre-cursors': 5,
      'number of pre-cursors in stock': 4,
      'average template occurrence': 909.8,
    },
    'metadata': {
      'created_at_iteration': 98,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>CC(=O)[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '9acc680c930290c191eb4a208e5cf64573c699d51b010ffd4ddde6d301764851',
          'classification': '0.0 Unrecognized',
          'library_occurence': 19,
          'policy_probability': 0.0004,
          'policy_probability_rank': 29,
          'policy_name': 'uspto',
          'template_code': 25786,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-C(=O)-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][C:27]([CH3:28])=[O:29]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'CC(C)=O',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': '01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 1107,
                  'policy_probability': 0.622600019,
                  'policy_probability_rank': 0,
                  'policy_name': 'uspto',
                  'template_code': 248,
                  'template': '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1.[OH:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[OH:21]',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'Cc1ccc(O)c(O)c1',
                    'is_chemical': true,
                    'in_stock': true,
                  },
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Cl)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 1863,
                          'policy_probability': 0.4481999874,
                          'policy_probability_rank': 0,
                          'policy_name': 'uspto',
                          'template_code': 4355,
                          'template': '[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:27])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:28])[cH:23][cH:24][c:25]2[cH:26]1.[Cl:27][S:29]([Cl:30])=[O:31]',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'O=S(Cl)Cl',
                            'is_chemical': true,
                            'in_stock': true,
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[O:1][C:2]([CH3:3])=[O:4]>>C[O:1][C:2]([CH3:3])=[O:4]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '110b771d77121b19f0a012d7c67263575768c8985321b84bf1629ea79c879d51',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 21925,
                                  'policy_probability': 0.4129999876,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 2818,
                                  'template': '[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:28])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COC(=O)C(C)c1ccc2cc(OC)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[C:1][C:6]([C:4]([O:3][CH3:2])=[O:5])[cH3:7]>>I[C:1].[CH3:2][O:3][C:4](=[O:5])[C:6][cH3:7]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '3c075bd2dba9541386d3ade4d2fb72e3c039b52e5142cf25a42d6cba6563011d',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 148,
                                          'policy_probability': 0.5683000088,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 10055,
                                          'template': '[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH;D3;+0:6](-[CH3;D1;+0:1])-[c:7]>>I-[CH3;D1;+0:1].[C;D1;H3:2]-[#8:3]-[C:4](=[O;D1;H0:5])-[CH2;D2;+0:6]-[c:7]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:10][I:13].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH2:9][C:11](=[O:12])[O:28][CH3:29])[cH:23][cH:24][c:25]2[cH:26]1',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'CI',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COC(=O)Cc1ccc2cc(OC)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7734470711,
      'number of reactions': 5,
      'number of pre-cursors': 5,
      'number of pre-cursors in stock': 4,
      'average template occurrence': 5012.4,
    },
    'metadata': {
      'created_at_iteration': 13,
      'is_solved': false,
    },
  },
  {
    'type': 'mol',
    'hide': false,
    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC)ccc2c1',
    'is_chemical': true,
    'in_stock': false,
    'children': [
      {
        'type': 'reaction',
        'hide': false,
        'smiles': '[C:1][O:2][cH3:3]>>CC(=O)[C:1].[O:2][cH3:3]',
        'is_reaction': true,
        'metadata': {
          'template_hash': '9acc680c930290c191eb4a208e5cf64573c699d51b010ffd4ddde6d301764851',
          'classification': '0.0 Unrecognized',
          'library_occurence': 19,
          'policy_probability': 0.0004,
          'policy_probability_rank': 29,
          'policy_name': 'uspto',
          'template_code': 25786,
          'template': '[CH3;D1;+0:1]-[O;H0;D2;+0:2]-[c:3]>>C-C(=O)-[CH3;D1;+0:1].[OH;D1;+0:2]-[c:3]',
          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][CH3:22])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1.[CH3:22][C:27]([CH3:28])=[O:29]',
        },
        'children': [
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'CC(C)=O',
            'is_chemical': true,
            'in_stock': true,
          },
          {
            'type': 'mol',
            'hide': false,
            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3O)ccc2c1',
            'is_chemical': true,
            'in_stock': false,
            'children': [
              {
                'type': 'reaction',
                'hide': false,
                'smiles': '[O:1][cH3:2]>>CC(=O)[O:1][cH3:2]',
                'is_reaction': true,
                'metadata': {
                  'template_hash': 'ba9cf922694c22786e94d6562607c5b880cd83d9f9ace2090b7b33f921f84f06',
                  'classification': '0.0 Unrecognized',
                  'library_occurence': 964,
                  'policy_probability': 0.0196000002,
                  'policy_probability_rank': 4,
                  'policy_name': 'uspto',
                  'template_code': 31122,
                  'template': '[OH;D1;+0:1]-[c:2]>>C-C(=O)-[O;H0;D2;+0:1]-[c:2]',
                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[OH:21])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1',
                },
                'children': [
                  {
                    'type': 'mol',
                    'hide': false,
                    'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3OC(C)=O)ccc2c1',
                    'is_chemical': true,
                    'in_stock': false,
                    'children': [
                      {
                        'type': 'reaction',
                        'hide': false,
                        'smiles': '[O:1]([C:3]([CH3:2])=[O:4])[c:5]([cH2:6])[cH2:7]>>O=C(O[O:1])c1cccc(Cl)c1.[CH3:2][C:3](=[O:4])[c:5]([cH2:6])[cH2:7]',
                        'is_reaction': true,
                        'metadata': {
                          'template_hash': '1ec0ea5b4cb26ec0598139be26b1ffece0af34c9422e1d0ac10d34892739d5c8',
                          'classification': '0.0 Unrecognized',
                          'library_occurence': 104,
                          'policy_probability': 0.0810000002,
                          'policy_probability_rank': 3,
                          'policy_name': 'uspto',
                          'template_code': 5183,
                          'template': '[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[O;H0;D2;+0:1]-[c;H0;D3;+0:5](:[c:6]):[c:7]>>Cl-c1:c:c:c:c(-C(=O)-O-[OH;D1;+0:1]):c:1.[C;D1;H3:2]-[C;H0;D3;+0:3](=[O;D1;H0:4])-[c;H0;D3;+0:5](:[c:6]):[c:7]',
                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[O:21][C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1.[OH:21][O:22][C:30]([c:31]1[cH:32][cH:33][cH:34][c:35]([Cl:36])[cH:37]1)=[O:38]',
                        },
                        'children': [
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'O=C(OO)c1cccc(Cl)c1',
                            'is_chemical': true,
                            'in_stock': true,
                          },
                          {
                            'type': 'mol',
                            'hide': false,
                            'smiles': 'COc1ccc2cc(C(C)C(=O)Oc3ccc(C)cc3C(C)=O)ccc2c1',
                            'is_chemical': true,
                            'in_stock': false,
                            'children': [
                              {
                                'type': 'reaction',
                                'hide': false,
                                'smiles': '[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]',
                                'is_reaction': true,
                                'metadata': {
                                  'template_hash': '01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26',
                                  'classification': '0.0 Unrecognized',
                                  'library_occurence': 1107,
                                  'policy_probability': 0.6410999894,
                                  'policy_probability_rank': 0,
                                  'policy_name': 'uspto',
                                  'template_code': 248,
                                  'template': '[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]',
                                  'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[O:13][c:14]3[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]3[C:27]([CH3:28])=[O:29])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:30])[cH:23][cH:24][c:25]2[cH:26]1.[OH:13][c:14]1[cH:15][cH:16][c:17]([CH3:18])[cH:19][c:20]1[C:27]([CH3:28])=[O:29]',
                                },
                                'children': [
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'CC(=O)c1cc(C)ccc1O',
                                    'is_chemical': true,
                                    'in_stock': true,
                                  },
                                  {
                                    'type': 'mol',
                                    'hide': false,
                                    'smiles': 'COc1ccc2cc(C(C)C(=O)Cl)ccc2c1',
                                    'is_chemical': true,
                                    'in_stock': false,
                                    'children': [
                                      {
                                        'type': 'reaction',
                                        'hide': false,
                                        'smiles': '[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]',
                                        'is_reaction': true,
                                        'metadata': {
                                          'template_hash': '19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9',
                                          'classification': '0.0 Unrecognized',
                                          'library_occurence': 1863,
                                          'policy_probability': 0.4481999874,
                                          'policy_probability_rank': 0,
                                          'policy_name': 'uspto',
                                          'template_code': 4355,
                                          'template': '[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]',
                                          'mapped_reaction_smiles': '[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[Cl:30])[cH:23][cH:24][c:25]2[cH:26]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[cH:7][c:8]([CH:9]([CH3:10])[C:11](=[O:12])[OH:27])[cH:23][cH:24][c:25]2[cH:26]1.[Cl:30][S:31]([Cl:32])=[O:33]',
                                        },
                                        'children': [
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'O=S(Cl)Cl',
                                            'is_chemical': true,
                                            'in_stock': true,
                                          },
                                          {
                                            'type': 'mol',
                                            'hide': false,
                                            'smiles': 'COc1ccc2cc(C(C)C(=O)O)ccc2c1',
                                            'is_chemical': true,
                                            'in_stock': false,
                                          },
                                        ],
                                      },
                                    ],
                                  },
                                ],
                              },
                            ],
                          },
                        ],
                      },
                    ],
                  },
                ],
              },
            ],
          },
        ],
      },
    ],
    'scores': {
      'state score': 0.7734470711,
      'number of reactions': 5,
      'number of pre-cursors': 5,
      'number of pre-cursors in stock': 4,
      'average template occurrence': 811.4,
    },
    'metadata': {
      'created_at_iteration': 99,
      'is_solved': false,
    },
  },
];


export const DEMO_DATA = [
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 1107,
                                  "policy_probability": 0.5583999753,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 248,
                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "O=C(Cl)C1CCCN1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1863,
                                                  "policy_probability": 0.0886000022,
                                                  "policy_probability_rank": 3,
                                                  "policy_name": "uspto",
                                                  "template_code": 4355,
                                                  "template": "[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]",
                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30]>>[Cl:30][S:31]([Cl:32])=[O:33].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(O)C1CCCN1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[CH3:1][C:2](=[O:3])[O:4]>>c1ccc(C[O:4][C:2]([CH3:1])=[O:3])cc1",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "2dfbe5ebafd5345589a283a7411e150fe2a301c7fc97f94643e2f41efc334c4e",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 2315,
                                                                  "policy_probability": 0.0286999997,
                                                                  "policy_probability_rank": 6,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 7727,
                                                                  "template": "[C:1]-[C:2](=[O;D1;H0:3])-[OH;D1;+0:4]>>[C:1]-[C:2](=[O;D1;H0:3])-[O;H0;D2;+0:4]-C-c1:c:c:c:c:c:1",
                                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "O=C(OCc1ccccc1)C1CCCN1",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[N:1]([CH3:2])[CH2:3][C:4](=[O:5])[O:6][CH2:7][cH3:8]>>CC(C)(C)OC(=O)[N:1]([CH3:2])[CH2:3][C:4](=[O:5])[O:6][CH2:7][cH3:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "9e740178f64e08553141c0f46ee6cbac3bc2bb929ad58b80bc07c0b96497234b",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 7,
                                                                                  "policy_probability": 0.3172999918,
                                                                                  "policy_probability_rank": 1,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 26415,
                                                                                  "template": "[C:2]-[NH;D2;+0:1]-[C:3]-[C:4](=[O;D1;H0:5])-[#8:6]-[C:7]-[c:8]>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]-[C:4](=[O;D1;H0:5])-[#8:6]-[C:7]-[c:8]",
                                                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1>>[N:13]1([C:28]([O:29][C:30]([CH3:31])([CH3:32])[CH3:33])=[O:34])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)OC(=O)N1CCCC1C(=O)OCc1ccccc1",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": false,
                                                                                      "children": [
                                                                                          {
                                                                                              "type": "reaction",
                                                                                              "hide": false,
                                                                                              "smiles": "[C:1]([cH3:2])[O:5][CH:4]=[O:3]>>Br[C:1][cH3:2].[O:3]=[CH:4][O:5]",
                                                                                              "is_reaction": true,
                                                                                              "metadata": {
                                                                                                  "template_hash": "6e0a0fa35ea74fe6ec010ab2dc1d62f9dc9e7138d23e03eeeefbda3de6c41be3",
                                                                                                  "classification": "0.0 Unrecognized",
                                                                                                  "library_occurence": 1212,
                                                                                                  "policy_probability": 0.4316000044,
                                                                                                  "policy_probability_rank": 0,
                                                                                                  "policy_name": "uspto",
                                                                                                  "template_code": 18357,
                                                                                                  "template": "[O;D1;H0:3]=[C:4]-[O;H0;D2;+0:5]-[CH2;D2;+0:1]-[c:2]>>Br-[CH2;D2;+0:1]-[c:2].[O;D1;H0:3]=[C:4]-[OH;D1;+0:5]",
                                                                                                  "mapped_reaction_smiles": "[N:13]1([C:28]([O:29][C:30]([CH3:31])([CH3:32])[CH3:33])=[O:34])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1>>[CH2:21]([c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1)[Br:35].[N:13]1([C:28]([O:29][C:30]([CH3:31])([CH3:32])[CH3:33])=[O:34])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                                                                              },
                                                                                              "children": [
                                                                                                  {
                                                                                                      "type": "mol",
                                                                                                      "hide": false,
                                                                                                      "smiles": "BrCc1ccccc1",
                                                                                                      "is_chemical": true,
                                                                                                      "in_stock": true
                                                                                                  },
                                                                                                  {
                                                                                                      "type": "mol",
                                                                                                      "hide": false,
                                                                                                      "smiles": "CC(C)(C)OC(=O)N1CCCC1C(=O)O",
                                                                                                      "is_chemical": true,
                                                                                                      "in_stock": false
                                                                                                  }
                                                                                              ]
                                                                                          }
                                                                                      ]
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=S(Cl)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[O:1][C:2]([CH3:3])=[O:4]>>CC(C)(C)[O:1][C:2]([CH3:3])=[O:4]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "b9a910a8443018b2725541937f75e6c441e220edbddaf92fde7cc76e7d6c0531",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 7512,
                                  "policy_probability": 0.0109999999,
                                  "policy_probability_rank": 4,
                                  "policy_name": "uspto",
                                  "template_code": 30956,
                                  "template": "[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-C(-C)(-C)-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc2c(c1)c(CC(=O)OC(C)(C)C)c(C)n2C(=O)c1ccc(Cl)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[n:5]([cH2:4])[cH2:6]>>Cl[C:1](=[O:2])[cH3:3].[cH2:4][n:5][cH2:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "6d581e03be0fd2c893db25c4afa6b3a9a79925f860fe3421d24c5531a882381a",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 246,
                                                  "policy_probability": 0.7189000249,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 18244,
                                                  "template": "[O;D1;H0:2]=[C;H0;D3;+0:1](-[c:3])-[n;H0;D3;+0:5](:[c:4]):[c:6]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[c:4]:[nH;D2;+0:5]:[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc2[nH]c(C)c(CC(=O)OC(C)(C)C)c2c1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[N:1][cH:2][c:3].O=[CH:4][C:5]",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "ca970095f1015377b40a1d2881e3a0fdbaf93b3ea66d7273385d9ac376eb583a",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 607,
                                                                  "policy_probability": 0.2741000056,
                                                                  "policy_probability_rank": 0,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 33696,
                                                                  "template": "[c;H0;D3;+0:3]1:[c;H0;D3;+0:2]:[nH;D2;+0:1]:[c;H0;D3;+0:4]:[c;H0;D3;+0:5]:1>>N-[NH;D2;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2>>[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:33])[cH:7][cH:8]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "COc1ccc(NN)cc1",
                                                                      "is_chemical": true,
                                                                      "in_stock": true
                                                                  },
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "CC(=O)CCC(=O)OC(C)(C)C",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:8][C:5]([CH3:4])([CH3:6])[CH3:7]>>O=[C:1]([CH3:2])[O:3].[CH3:4][C:5]([CH3:6])([CH3:7])[O:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "cfe38639da8fcfe3ed728df8b13d171dd8a4279bcb6a0839ebe0f753e7996b95",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 189,
                                                                                  "policy_probability": 0.1999000013,
                                                                                  "policy_probability_rank": 0,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 34584,
                                                                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;H0;D1;+0:3])-[O;H0;D2;+0:8]-[C:5](-[C;D1;H3:4])(-[C;D1;H3:6])-[C;D1;H3:7]>>O=[C;H0;D3;+0:1](-[C:2])-[OH;D1;+0:3].[C;D1;H3:4]-[C:5](-[C;D1;H3:6])(-[C;D1;H3:7])-[OH;D1;+0:8]",
                                                                                  "mapped_reaction_smiles": "[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47]>>[CH2:9]([CH2:10][C:11]([OH:12])=[O:48])[C:30]([CH3:31])=[O:47].[OH:42][C:43]([CH3:44])([CH3:45])[CH3:46]"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  },
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(=O)CCC(=O)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8372101461,
          "number of reactions": 10,
          "number of pre-cursors": 8,
          "number of pre-cursors in stock": 7,
          "average template occurrence": 2873.1
      },
      "metadata": {
          "created_at_iteration": 79,
          "is_solved": false
      }
  },
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 1107,
                                  "policy_probability": 0.5583999753,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 248,
                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "O=C(Cl)C1CCCN1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1863,
                                                  "policy_probability": 0.0886000022,
                                                  "policy_probability_rank": 3,
                                                  "policy_name": "uspto",
                                                  "template_code": 4355,
                                                  "template": "[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]",
                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30]>>[Cl:30][S:31]([Cl:32])=[O:33].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(O)C1CCCN1",
                                                      "is_chemical": true,
                                                      "in_stock": false
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=S(Cl)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[NH:1][cH:2][c:3].O=[CH:4][C:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "ec2c35151075985c7538d160ebafb7899a46e50cfeae3cf3a8bc0677a4cc1f07",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 115,
                                  "policy_probability": 0.7286000252,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 39261,
                                  "template": "[c;H0;D3;+0:2]1:[c;H0;D3;+0:3]:[c;H0;D3;+0:5]:[c;H0;D3;+0:4]:[n;H0;D3;+0:1]:1>>N-[N;H0;D3;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH2:9]([CH2:10][C:11](=[O:12])[OH:42])[C:30]([CH3:31])=[O:43].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([N:32]([C:33](=[O:34])[c:35]2[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]2)[NH2:44])[cH:7][cH:8]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "CC(=O)CCC(=O)O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc(N(N)C(=O)c2ccc(Cl)cc2)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[N:5]([NH2:4])[cH3:6]>>Cl[C:1](=[O:2])[cH3:3].[NH2:4][N:5][cH3:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "4305b703934eb792fa3f012faa9fb14dcd6f579babef07e115ee93c013311166",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 10,
                                                  "policy_probability": 0.6204000115,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 11243,
                                                  "template": "[N;D1;H2:4]-[N;H0;D3;+0:5](-[c:6])-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[N;D1;H2:4]-[NH;D2;+0:5]-[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([N:32]([C:33](=[O:34])[c:35]2[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]2)[NH2:44])[cH:7][cH:8]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:42].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:44])[cH:7][cH:8]1"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc(NN)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8282195956,
          "number of reactions": 5,
          "number of pre-cursors": 6,
          "number of pre-cursors in stock": 5,
          "average template occurrence": 3353.6
      },
      "metadata": {
          "created_at_iteration": 3,
          "is_solved": false
      }
  },
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 1107,
                                  "policy_probability": 0.5583999753,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 248,
                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "O=C(Cl)C1CCCN1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1863,
                                                  "policy_probability": 0.0886000022,
                                                  "policy_probability_rank": 3,
                                                  "policy_name": "uspto",
                                                  "template_code": 4355,
                                                  "template": "[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]",
                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30]>>[Cl:30][S:31]([Cl:32])=[O:33].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(O)C1CCCN1",
                                                      "is_chemical": true,
                                                      "in_stock": false
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=S(Cl)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[O:1][C:2]([CH3:3])=[O:4]>>CC(C)(C)[O:1][C:2]([CH3:3])=[O:4]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "b9a910a8443018b2725541937f75e6c441e220edbddaf92fde7cc76e7d6c0531",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 7512,
                                  "policy_probability": 0.0109999999,
                                  "policy_probability_rank": 4,
                                  "policy_name": "uspto",
                                  "template_code": 30956,
                                  "template": "[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-C(-C)(-C)-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc2c(c1)c(CC(=O)OC(C)(C)C)c(C)n2C(=O)c1ccc(Cl)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[n:5]([cH2:4])[cH2:6]>>Cl[C:1](=[O:2])[cH3:3].[cH2:4][n:5][cH2:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "6d581e03be0fd2c893db25c4afa6b3a9a79925f860fe3421d24c5531a882381a",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 246,
                                                  "policy_probability": 0.7189000249,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 18244,
                                                  "template": "[O;D1;H0:2]=[C;H0;D3;+0:1](-[c:3])-[n;H0;D3;+0:5](:[c:4]):[c:6]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[c:4]:[nH;D2;+0:5]:[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc2[nH]c(C)c(CC(=O)OC(C)(C)C)c2c1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[N:1][cH:2][c:3].O=[CH:4][C:5]",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "ca970095f1015377b40a1d2881e3a0fdbaf93b3ea66d7273385d9ac376eb583a",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 607,
                                                                  "policy_probability": 0.2741000056,
                                                                  "policy_probability_rank": 0,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 33696,
                                                                  "template": "[c;H0;D3;+0:3]1:[c;H0;D3;+0:2]:[nH;D2;+0:1]:[c;H0;D3;+0:4]:[c;H0;D3;+0:5]:1>>N-[NH;D2;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2>>[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:33])[cH:7][cH:8]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "COc1ccc(NN)cc1",
                                                                      "is_chemical": true,
                                                                      "in_stock": true
                                                                  },
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "CC(=O)CCC(=O)OC(C)(C)C",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:8][C:5]([CH3:4])([CH3:6])[CH3:7]>>O=[C:1]([CH3:2])[O:3].[CH3:4][C:5]([CH3:6])([CH3:7])[O:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "cfe38639da8fcfe3ed728df8b13d171dd8a4279bcb6a0839ebe0f753e7996b95",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 189,
                                                                                  "policy_probability": 0.1999000013,
                                                                                  "policy_probability_rank": 0,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 34584,
                                                                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;H0;D1;+0:3])-[O;H0;D2;+0:8]-[C:5](-[C;D1;H3:4])(-[C;D1;H3:6])-[C;D1;H3:7]>>O=[C;H0;D3;+0:1](-[C:2])-[OH;D1;+0:3].[C;D1;H3:4]-[C:5](-[C;D1;H3:6])(-[C;D1;H3:7])-[OH;D1;+0:8]",
                                                                                  "mapped_reaction_smiles": "[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47]>>[CH2:9]([CH2:10][C:11]([OH:12])=[O:48])[C:30]([CH3:31])=[O:47].[OH:42][C:43]([CH3:44])([CH3:45])[CH3:46]"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  },
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(=O)CCC(=O)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8277327854,
          "number of reactions": 7,
          "number of pre-cursors": 7,
          "number of pre-cursors in stock": 6,
          "average template occurrence": 3599.5714285714
      },
      "metadata": {
          "created_at_iteration": 79,
          "is_solved": false
      }
  },
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 1107,
                                  "policy_probability": 0.5583999753,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 248,
                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "O=C(Cl)C1CCCN1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1863,
                                                  "policy_probability": 0.0886000022,
                                                  "policy_probability_rank": 3,
                                                  "policy_name": "uspto",
                                                  "template_code": 4355,
                                                  "template": "[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]",
                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30]>>[Cl:30][S:31]([Cl:32])=[O:33].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(O)C1CCCN1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[CH3:1][C:2](=[O:3])[O:4]>>c1ccc(C[O:4][C:2]([CH3:1])=[O:3])cc1",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "2dfbe5ebafd5345589a283a7411e150fe2a301c7fc97f94643e2f41efc334c4e",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 2315,
                                                                  "policy_probability": 0.0286999997,
                                                                  "policy_probability_rank": 6,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 7727,
                                                                  "template": "[C:1]-[C:2](=[O;D1;H0:3])-[OH;D1;+0:4]>>[C:1]-[C:2](=[O;D1;H0:3])-[O;H0;D2;+0:4]-C-c1:c:c:c:c:c:1",
                                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "O=C(OCc1ccccc1)C1CCCN1",
                                                                      "is_chemical": true,
                                                                      "in_stock": false
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=S(Cl)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[O:1][C:2]([CH3:3])=[O:4]>>CC(C)(C)[O:1][C:2]([CH3:3])=[O:4]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "b9a910a8443018b2725541937f75e6c441e220edbddaf92fde7cc76e7d6c0531",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 7512,
                                  "policy_probability": 0.0109999999,
                                  "policy_probability_rank": 4,
                                  "policy_name": "uspto",
                                  "template_code": 30956,
                                  "template": "[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-C(-C)(-C)-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc2c(c1)c(CC(=O)OC(C)(C)C)c(C)n2C(=O)c1ccc(Cl)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[n:5]([cH2:4])[cH2:6]>>Cl[C:1](=[O:2])[cH3:3].[cH2:4][n:5][cH2:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "6d581e03be0fd2c893db25c4afa6b3a9a79925f860fe3421d24c5531a882381a",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 246,
                                                  "policy_probability": 0.7189000249,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 18244,
                                                  "template": "[O;D1;H0:2]=[C;H0;D3;+0:1](-[c:3])-[n;H0;D3;+0:5](:[c:4]):[c:6]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[c:4]:[nH;D2;+0:5]:[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc2[nH]c(C)c(CC(=O)OC(C)(C)C)c2c1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[N:1][cH:2][c:3].O=[CH:4][C:5]",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "ca970095f1015377b40a1d2881e3a0fdbaf93b3ea66d7273385d9ac376eb583a",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 607,
                                                                  "policy_probability": 0.2741000056,
                                                                  "policy_probability_rank": 0,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 33696,
                                                                  "template": "[c;H0;D3;+0:3]1:[c;H0;D3;+0:2]:[nH;D2;+0:1]:[c;H0;D3;+0:4]:[c;H0;D3;+0:5]:1>>N-[NH;D2;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2>>[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:33])[cH:7][cH:8]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "COc1ccc(NN)cc1",
                                                                      "is_chemical": true,
                                                                      "in_stock": true
                                                                  },
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "CC(=O)CCC(=O)OC(C)(C)C",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:8][C:5]([CH3:4])([CH3:6])[CH3:7]>>O=[C:1]([CH3:2])[O:3].[CH3:4][C:5]([CH3:6])([CH3:7])[O:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "cfe38639da8fcfe3ed728df8b13d171dd8a4279bcb6a0839ebe0f753e7996b95",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 189,
                                                                                  "policy_probability": 0.1999000013,
                                                                                  "policy_probability_rank": 0,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 34584,
                                                                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;H0;D1;+0:3])-[O;H0;D2;+0:8]-[C:5](-[C;D1;H3:4])(-[C;D1;H3:6])-[C;D1;H3:7]>>O=[C;H0;D3;+0:1](-[C:2])-[OH;D1;+0:3].[C;D1;H3:4]-[C:5](-[C;D1;H3:6])(-[C;D1;H3:7])-[OH;D1;+0:8]",
                                                                                  "mapped_reaction_smiles": "[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47]>>[CH2:9]([CH2:10][C:11]([OH:12])=[O:48])[C:30]([CH3:31])=[O:47].[OH:42][C:43]([CH3:44])([CH3:45])[CH3:46]"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  },
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(=O)CCC(=O)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8277327854,
          "number of reactions": 8,
          "number of pre-cursors": 7,
          "number of pre-cursors in stock": 6,
          "average template occurrence": 3439.0
      },
      "metadata": {
          "created_at_iteration": 79,
          "is_solved": false
      }
  },
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 1107,
                                  "policy_probability": 0.5583999753,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 248,
                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "O=C(Cl)C1CCCN1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[Cl:1][C:2]([CH3:3])=[O:4]>>O=S(Cl)[Cl:1].O[C:2]([CH3:3])=[O:4]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "19f2a0f4c6a46a956748171b40376b6f742843423bcde2abbd8b3a6cf6b0a7b9",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1863,
                                                  "policy_probability": 0.0886000022,
                                                  "policy_probability_rank": 3,
                                                  "policy_name": "uspto",
                                                  "template_code": 4355,
                                                  "template": "[C:3]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:4]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](-[C:3])=[O;D1;H0:4]",
                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:30]>>[Cl:30][S:31]([Cl:32])=[O:33].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(O)C1CCCN1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[CH3:1][C:2](=[O:3])[O:4]>>c1ccc(C[O:4][C:2]([CH3:1])=[O:3])cc1",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "2dfbe5ebafd5345589a283a7411e150fe2a301c7fc97f94643e2f41efc334c4e",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 2315,
                                                                  "policy_probability": 0.0286999997,
                                                                  "policy_probability_rank": 6,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 7727,
                                                                  "template": "[C:1]-[C:2](=[O;D1;H0:3])-[OH;D1;+0:4]>>[C:1]-[C:2](=[O;D1;H0:3])-[O;H0;D2;+0:4]-C-c1:c:c:c:c:c:1",
                                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:20]>>[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "O=C(OCc1ccccc1)C1CCCN1",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[N:1]([CH3:2])[CH2:3][C:4](=[O:5])[O:6][CH2:7][cH3:8]>>CC(C)(C)OC(=O)[N:1]([CH3:2])[CH2:3][C:4](=[O:5])[O:6][CH2:7][cH3:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "9e740178f64e08553141c0f46ee6cbac3bc2bb929ad58b80bc07c0b96497234b",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 7,
                                                                                  "policy_probability": 0.3172999918,
                                                                                  "policy_probability_rank": 1,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 26415,
                                                                                  "template": "[C:2]-[NH;D2;+0:1]-[C:3]-[C:4](=[O;D1;H0:5])-[#8:6]-[C:7]-[c:8]>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]-[C:4](=[O;D1;H0:5])-[#8:6]-[C:7]-[c:8]",
                                                                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1>>[N:13]1([C:28]([O:29][C:30]([CH3:31])([CH3:32])[CH3:33])=[O:34])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][CH2:21][c:22]1[cH:23][cH:24][cH:25][cH:26][cH:27]1"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)OC(=O)N1CCCC1C(=O)OCc1ccccc1",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": false
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=S(Cl)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[O:1][C:2]([CH3:3])=[O:4]>>CC(C)(C)[O:1][C:2]([CH3:3])=[O:4]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "b9a910a8443018b2725541937f75e6c441e220edbddaf92fde7cc76e7d6c0531",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 7512,
                                  "policy_probability": 0.0109999999,
                                  "policy_probability_rank": 4,
                                  "policy_name": "uspto",
                                  "template_code": 30956,
                                  "template": "[C:3]-[C:2](=[O;D1;H0:4])-[OH;D1;+0:1]>>C-C(-C)(-C)-[O;H0;D2;+0:1]-[C:2](-[C:3])=[O;D1;H0:4]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc2c(c1)c(CC(=O)OC(C)(C)C)c(C)n2C(=O)c1ccc(Cl)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[n:5]([cH2:4])[cH2:6]>>Cl[C:1](=[O:2])[cH3:3].[cH2:4][n:5][cH2:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "6d581e03be0fd2c893db25c4afa6b3a9a79925f860fe3421d24c5531a882381a",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 246,
                                                  "policy_probability": 0.7189000249,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 18244,
                                                  "template": "[O;D1;H0:2]=[C;H0;D3;+0:1](-[c:3])-[n;H0;D3;+0:5](:[c:4]):[c:6]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[c:4]:[nH;D2;+0:5]:[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc2[nH]c(C)c(CC(=O)OC(C)(C)C)c2c1",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[N:1][cH:2][c:3].O=[CH:4][C:5]",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "ca970095f1015377b40a1d2881e3a0fdbaf93b3ea66d7273385d9ac376eb583a",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 607,
                                                                  "policy_probability": 0.2741000056,
                                                                  "policy_probability_rank": 0,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 33696,
                                                                  "template": "[c;H0;D3;+0:3]1:[c;H0;D3;+0:2]:[nH;D2;+0:1]:[c;H0;D3;+0:4]:[c;H0;D3;+0:5]:1>>N-[NH;D2;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[c:30]([CH3:31])[nH:32]2>>[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:33])[cH:7][cH:8]1"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "COc1ccc(NN)cc1",
                                                                      "is_chemical": true,
                                                                      "in_stock": true
                                                                  },
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "CC(=O)CCC(=O)OC(C)(C)C",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:8][C:5]([CH3:4])([CH3:6])[CH3:7]>>O=[C:1]([CH3:2])[O:3].[CH3:4][C:5]([CH3:6])([CH3:7])[O:8]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "cfe38639da8fcfe3ed728df8b13d171dd8a4279bcb6a0839ebe0f753e7996b95",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 189,
                                                                                  "policy_probability": 0.1999000013,
                                                                                  "policy_probability_rank": 0,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 34584,
                                                                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;H0;D1;+0:3])-[O;H0;D2;+0:8]-[C:5](-[C;D1;H3:4])(-[C;D1;H3:6])-[C;D1;H3:7]>>O=[C;H0;D3;+0:1](-[C:2])-[OH;D1;+0:3].[C;D1;H3:4]-[C:5](-[C;D1;H3:6])(-[C;D1;H3:7])-[OH;D1;+0:8]",
                                                                                  "mapped_reaction_smiles": "[CH2:9]([CH2:10][C:11](=[O:12])[O:42][C:43]([CH3:44])([CH3:45])[CH3:46])[C:30]([CH3:31])=[O:47]>>[CH2:9]([CH2:10][C:11]([OH:12])=[O:48])[C:30]([CH3:31])=[O:47].[OH:42][C:43]([CH3:44])([CH3:45])[CH3:46]"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  },
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(=O)CCC(=O)O",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8277327854,
          "number of reactions": 9,
          "number of pre-cursors": 7,
          "number of pre-cursors in stock": 6,
          "average template occurrence": 3057.6666666667
      },
      "metadata": {
          "created_at_iteration": 79,
          "is_solved": false
      }
  },
  {
      "type": "mol",
      "hide": false,
      "smiles": "COc1ccc2c(c1)c(CC(=O)N1CCCC1C(=O)Oc1ccc(C)cc1OC)c(C)n2C(=O)c1ccc(Cl)cc1",
      "is_chemical": true,
      "in_stock": false,
      "children": [
          {
              "type": "reaction",
              "hide": false,
              "smiles": "[C:1]([CH3:2])(=[O:3])[N:5]([CH3:4])[CH3:6]>>O[C:1]([CH3:2])=[O:3].[CH3:4][N:5][CH3:6]",
              "is_reaction": true,
              "metadata": {
                  "template_hash": "46395a26b4b93170d8898450d8ec58e43080e87af1ad5c1330db3a944d0d57ac",
                  "classification": "0.0 Unrecognized",
                  "library_occurence": 13673,
                  "policy_probability": 0.0048000002,
                  "policy_probability_rank": 4,
                  "policy_name": "uspto",
                  "template_code": 11766,
                  "template": "[C:4]-[N;H0;D3;+0:5](-[C:6])-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3]>>O-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[C:4]-[NH;D2;+0:5]-[C:6]",
                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[N:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1.[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
              },
              "children": [
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[N:1]([CH3:2])[CH2:3][CH:4]=[O:5]>>CC(C)(C)OC(=O)[N:1]([CH3:2])[CH2:3][CH:4]=[O:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "455b92910bd75e5ccc17456b6401025eaa099933fec292c5b20cb3846588d20c",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 49,
                                  "policy_probability": 0.0170000009,
                                  "policy_probability_rank": 6,
                                  "policy_name": "uspto",
                                  "template_code": 11627,
                                  "template": "[C:2]-[NH;D2;+0:1]-[C:3]-[C:4]=[O;D1;H0:5]>>C-C(-C)(-C)-O-C(=O)-[N;H0;D3;+0:1](-[C:2])-[C:3]-[C:4]=[O;D1;H0:5]",
                                  "mapped_reaction_smiles": "[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1cc(C)ccc1OC(=O)C1CCCN1C(=O)OC(C)(C)C",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1]([CH3:2])(=[O:3])[O:4][cH3:5]>>Cl[C:1]([CH3:2])=[O:3].[O:4][cH3:5]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "01643639d6a55c16f7f30c6505aeea5e206f45f41edb94fb92d12bf48278cc26",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 1107,
                                                  "policy_probability": 0.2809000015,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 248,
                                                  "template": "[C:2]-[C;H0;D3;+0:1](=[O;D1;H0:3])-[O;H0;D2;+0:4]-[c:5]>>Cl-[C;H0;D3;+0:1](-[C:2])=[O;D1;H0:3].[OH;D1;+0:4]-[c:5]",
                                                  "mapped_reaction_smiles": "[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[O:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]>>[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:37].[OH:20][c:21]1[cH:22][cH:23][c:24]([CH3:25])[cH:26][c:27]1[O:28][CH3:29]"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1cc(C)ccc1O",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "CC(C)(C)OC(=O)N1CCCC1C(=O)Cl",
                                                      "is_chemical": true,
                                                      "in_stock": false,
                                                      "children": [
                                                          {
                                                              "type": "reaction",
                                                              "hide": false,
                                                              "smiles": "[Cl:1][C:2](=[O:3])[CH2:4][NH2:5]>>O=S(Cl)[Cl:1].O[C:2](=[O:3])[CH2:4][NH2:5]",
                                                              "is_reaction": true,
                                                              "metadata": {
                                                                  "template_hash": "56bc62f38aa912192a49287110f90c95c2e62b13c7d795cf13cef46ee7c36b81",
                                                                  "classification": "0.0 Unrecognized",
                                                                  "library_occurence": 20,
                                                                  "policy_probability": 0.1403000057,
                                                                  "policy_probability_rank": 1,
                                                                  "policy_name": "uspto",
                                                                  "template_code": 14474,
                                                                  "template": "[#7:5]-[C:4]-[C;H0;D3;+0:2](-[Cl;H0;D1;+0:1])=[O;D1;H0:3]>>Cl-S(=O)-[Cl;H0;D1;+0:1].O-[C;H0;D3;+0:2](=[O;D1;H0:3])-[C:4]-[#7:5]",
                                                                  "mapped_reaction_smiles": "[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[Cl:37]>>[Cl:37][S:38]([Cl:39])=[O:40].[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:41]"
                                                              },
                                                              "children": [
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "O=S(Cl)Cl",
                                                                      "is_chemical": true,
                                                                      "in_stock": true
                                                                  },
                                                                  {
                                                                      "type": "mol",
                                                                      "hide": false,
                                                                      "smiles": "CC(C)(C)OC(=O)N1CCCC1C(=O)O",
                                                                      "is_chemical": true,
                                                                      "in_stock": false,
                                                                      "children": [
                                                                          {
                                                                              "type": "reaction",
                                                                              "hide": false,
                                                                              "smiles": "[C:1](=[O:2])([O:3][C:4]([CH3:5])([CH3:6])[CH3:7])[N:9]([CH3:8])[CH2:10][C:11](=[O:12])[OH:13]>>CC(C)(C)OC(=O)O[C:1](=[O:2])[O:3][C:4]([CH3:5])([CH3:6])[CH3:7].[CH3:8][N:9][CH2:10][C:11](=[O:12])[OH:13]",
                                                                              "is_reaction": true,
                                                                              "metadata": {
                                                                                  "template_hash": "20032004b7af43577f93b5d9f6b6511fa3c9c8a14315d4c007cebb48f1107c66",
                                                                                  "classification": "0.0 Unrecognized",
                                                                                  "library_occurence": 447,
                                                                                  "policy_probability": 0.695299983,
                                                                                  "policy_probability_rank": 0,
                                                                                  "policy_name": "uspto",
                                                                                  "template_code": 5386,
                                                                                  "template": "[C:8]-[N;H0;D3;+0:9](-[C:10]-[C:11](=[O;D1;H0:12])-[O;D1;H1:13])-[C;H0;D3;+0:1](=[O;D1;H0:2])-[#8:3]-[C:4](-[C;D1;H3:5])(-[C;D1;H3:6])-[C;D1;H3:7]>>C-C(-C)(-C)-O-C(=O)-O-[C;H0;D3;+0:1](=[O;D1;H0:2])-[#8:3]-[C:4](-[C;D1;H3:5])(-[C;D1;H3:6])-[C;D1;H3:7].[C:8]-[NH;D2;+0:9]-[C:10]-[C:11](=[O;D1;H0:12])-[O;D1;H1:13]",
                                                                                  "mapped_reaction_smiles": "[N:13]1([C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])=[O:36])[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:41]>>[C:30]([O:31][C:32]([CH3:33])([CH3:34])[CH3:35])(=[O:36])[O:37][C:38]([O:39][C:40]([CH3:42])([CH3:43])[CH3:44])=[O:45].[NH:13]1[CH2:14][CH2:15][CH2:16][CH:17]1[C:18](=[O:19])[OH:41]"
                                                                              },
                                                                              "children": [
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "O=C(O)C1CCCN1",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": false
                                                                                  },
                                                                                  {
                                                                                      "type": "mol",
                                                                                      "hide": false,
                                                                                      "smiles": "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C",
                                                                                      "is_chemical": true,
                                                                                      "in_stock": true
                                                                                  }
                                                                              ]
                                                                          }
                                                                      ]
                                                                  }
                                                              ]
                                                          }
                                                      ]
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  },
                  {
                      "type": "mol",
                      "hide": false,
                      "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                      "is_chemical": true,
                      "in_stock": false,
                      "children": [
                          {
                              "type": "reaction",
                              "hide": false,
                              "smiles": "[n:1]1[cH:2][cH:3][cH:5][cH:4]1>>N[NH:1][cH:2][c:3].O=[CH:4][C:5]",
                              "is_reaction": true,
                              "metadata": {
                                  "template_hash": "ec2c35151075985c7538d160ebafb7899a46e50cfeae3cf3a8bc0677a4cc1f07",
                                  "classification": "0.0 Unrecognized",
                                  "library_occurence": 115,
                                  "policy_probability": 0.7286000252,
                                  "policy_probability_rank": 0,
                                  "policy_name": "uspto",
                                  "template_code": 39261,
                                  "template": "[c;H0;D3;+0:2]1:[c;H0;D3;+0:3]:[c;H0;D3;+0:5]:[c;H0;D3;+0:4]:[n;H0;D3;+0:1]:1>>N-[N;H0;D3;+0:1]-[c;H0;D3;+0:2]:[cH;D2;+0:3].O=[C;H0;D3;+0:4]-[CH2;D2;+0:5]",
                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]2[c:7]([cH:8]1)[c:9]([CH2:10][C:11](=[O:12])[OH:42])[c:30]([CH3:31])[n:32]2[C:33](=[O:34])[c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1>>[CH2:9]([CH2:10][C:11](=[O:12])[OH:42])[C:30]([CH3:31])=[O:43].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([N:32]([C:33](=[O:34])[c:35]2[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]2)[NH2:44])[cH:7][cH:8]1"
                              },
                              "children": [
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "CC(=O)CCC(=O)O",
                                      "is_chemical": true,
                                      "in_stock": true
                                  },
                                  {
                                      "type": "mol",
                                      "hide": false,
                                      "smiles": "COc1ccc(N(N)C(=O)c2ccc(Cl)cc2)cc1",
                                      "is_chemical": true,
                                      "in_stock": false,
                                      "children": [
                                          {
                                              "type": "reaction",
                                              "hide": false,
                                              "smiles": "[C:1](=[O:2])([cH3:3])[N:5]([NH2:4])[cH3:6]>>Cl[C:1](=[O:2])[cH3:3].[NH2:4][N:5][cH3:6]",
                                              "is_reaction": true,
                                              "metadata": {
                                                  "template_hash": "4305b703934eb792fa3f012faa9fb14dcd6f579babef07e115ee93c013311166",
                                                  "classification": "0.0 Unrecognized",
                                                  "library_occurence": 10,
                                                  "policy_probability": 0.6204000115,
                                                  "policy_probability_rank": 0,
                                                  "policy_name": "uspto",
                                                  "template_code": 11243,
                                                  "template": "[N;D1;H2:4]-[N;H0;D3;+0:5](-[c:6])-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3]>>Cl-[C;H0;D3;+0:1](=[O;D1;H0:2])-[c:3].[N;D1;H2:4]-[NH;D2;+0:5]-[c:6]",
                                                  "mapped_reaction_smiles": "[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([N:32]([C:33](=[O:34])[c:35]2[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]2)[NH2:44])[cH:7][cH:8]1>>[C:33](=[O:34])([c:35]1[cH:36][cH:37][c:38]([Cl:39])[cH:40][cH:41]1)[Cl:42].[CH3:1][O:2][c:3]1[cH:4][cH:5][c:6]([NH:32][NH2:44])[cH:7][cH:8]1"
                                              },
                                              "children": [
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "COc1ccc(NN)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  },
                                                  {
                                                      "type": "mol",
                                                      "hide": false,
                                                      "smiles": "O=C(Cl)c1ccc(Cl)cc1",
                                                      "is_chemical": true,
                                                      "in_stock": true
                                                  }
                                              ]
                                          }
                                      ]
                                  }
                              ]
                          }
                      ]
                  }
              ]
          }
      ],
      "scores": {
          "state score": 0.8277327854,
          "number of reactions": 7,
          "number of pre-cursors": 7,
          "number of pre-cursors in stock": 6,
          "average template occurrence": 2203.0
      },
      "metadata": {
          "created_at_iteration": 46,
          "is_solved": false
      }
  }
]