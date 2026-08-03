/* External reference fixtures for Control comparisons — maintained by hand.
 *
 * Produced once with (numpy 2.4.3 | scipy 1.17.1 | statsmodels 0.14.6 | pandas 2.3.3, python 3.13):
 *   - scipy.stats.dunnett (seeded QMC)         -> Dunnett statistic + adjusted p
 *   - scipy.stats.ttest_ind(equal_var=False)   -> Welch t, df, raw p, CI on the mean diff
 *   - statsmodels multipletests(method='holm') -> Holm over the k-1 control-vs-each p's
 *
 * FIXTURES       - synthetic / published cases (control + treated arrays).
 * DEMOG_FIXTURES - real data; the test reads DEMOG_FILE (the platform's demo demog.csv) and runs
 *                  the named (group, feature, control, alpha) combos. Expected values are keyed by
 *                  group label. Only no-null group columns are used (RACE/SEX/DIS_POP); null
 *                  features (HEIGHT/WEIGHT/AGE) are dropped identically by factorize and pandas.
 *
 * To add cases, compute the references with the libraries above and paste them here.
 */

/* eslint-disable */

export const ALPHA = 0.05;
export const DEMOG_FILE = "System:DemoFiles/demog.csv";

export interface DunnettExpected {statistic: number; pValueAdj: number;}
export interface HolmWelchExpected {
  meanDiff: number; statistic: number; df: number;
  pValueRaw: number; pValueAdj: number; ciLow: number; ciHigh: number; hedgesG: number;
}
export interface ControlComparisonsFixture {
  name: string;
  control: number[];
  treated: number[][];
  dunnett: DunnettExpected[];
  holmWelch: HolmWelchExpected[];
}

export interface DunnettDemogExpected extends DunnettExpected {groupLabel: string;}
export interface HolmWelchDemogExpected extends HolmWelchExpected {groupLabel: string;}
export interface ControlComparisonsDemogFixture {
  name: string;
  groupCol: string;
  featureCol: string;
  control: string;
  alpha: number;
  dunnett: DunnettDemogExpected[];
  holmWelch: HolmWelchDemogExpected[];
}

export const FIXTURES: ControlComparisonsFixture[] = [
  {
    "name": "[Dunnett 1955] blood counts",
    "control": [
      7.4,
      8.5,
      7.2,
      8.24,
      9.84,
      8.32
    ],
    "treated": [
      [
        9.76,
        8.8,
        7.68,
        9.36
      ],
      [
        12.8,
        9.68,
        12.16,
        9.2,
        10.55
      ]
    ],
    "dunnett": [
      {
        "statistic": 0.8570321259343803,
        "pValueAdj": 0.6200702963110303
      },
      {
        "statistic": 3.693752311235726,
        "pValueAdj": 0.005764730919475047
      }
    ],
    "holmWelch": [
      {
        "meanDiff": 0.6499999999999986,
        "statistic": 1.0963745079232525,
        "df": 6.7758809316942905,
        "pValueRaw": 0.31034466666098365,
        "ciLow": -0.7613533565627968,
        "ciHigh": 2.061353356562794,
        "hedgesG": 0.6233611465235627,
        "pValueAdj": 0.31034466666098365
      },
      {
        "meanDiff": 2.6279999999999983,
        "statistic": 3.3053143622690477,
        "df": 6.333034350606968,
        "pValueRaw": 0.01505982907786063,
        "ciLow": 0.7070340249846427,
        "ciHigh": 4.548965975015354,
        "hedgesG": 1.7897225709422269,
        "pValueAdj": 0.03011965815572126
      }
    ]
  },
  {
    "name": "balanced, 3 treated",
    "control": [
      1.0,
      2.0,
      3.0,
      2.5,
      1.5
    ],
    "treated": [
      [
        3.0,
        4.0,
        5.0,
        4.5
      ],
      [
        5.0,
        6.0,
        7.0,
        6.5
      ],
      [
        7.0,
        8.0,
        9.0,
        8.5
      ]
    ],
    "dunnett": [
      {
        "statistic": 3.7940263593345884,
        "pValueAdj": 0.006049469035476207
      },
      {
        "statistic": 7.364874697531849,
        "pValueAdj": 8.057462626887002e-06
      },
      {
        "statistic": 10.93572303572911,
        "pValueAdj": 1.1575785319184462e-07
      }
    ],
    "holmWelch": [
      {
        "meanDiff": 2.125,
        "statistic": 3.833394742814427,
        "df": 6.302353651176825,
        "pValueRaw": 0.007852719505041753,
        "ciLow": 0.7842004551727213,
        "ciHigh": 3.465799544827279,
        "hedgesG": 2.260273609589547,
        "pValueAdj": 0.007852719505041753
      },
      {
        "meanDiff": 4.125,
        "statistic": 7.4412956772280054,
        "df": 6.302353651176825,
        "pValueRaw": 0.0002404659202038211,
        "ciLow": 2.784200455172721,
        "ciHigh": 5.465799544827279,
        "hedgesG": 4.387589948026768,
        "pValueAdj": 0.0004809318404076422
      },
      {
        "meanDiff": 6.125,
        "statistic": 11.049196611641584,
        "df": 6.302353651176825,
        "pValueRaw": 2.3353045069891925e-05,
        "ciLow": 4.784200455172721,
        "ciHigh": 7.465799544827279,
        "hedgesG": 6.514906286463989,
        "pValueAdj": 7.005913520967577e-05
      }
    ]
  },
  {
    "name": "unbalanced, 2 treated",
    "control": [
      10.0,
      11.0,
      9.5,
      10.5,
      10.2,
      11.3,
      9.8
    ],
    "treated": [
      [
        12.0,
        11.5,
        13.0
      ],
      [
        14.0,
        15.0,
        14.5,
        13.5
      ]
    ],
    "dunnett": [
      {
        "statistic": 3.9784013687041444,
        "pValueAdj": 0.0040349520872471345
      },
      {
        "statistic": 9.34453984503784,
        "pValueAdj": 8.659753906181678e-07
      }
    ],
    "holmWelch": [
      {
        "meanDiff": 1.8380952380952387,
        "statistic": 3.6449171884618377,
        "df": 3.3162790406416565,
        "pValueRaw": 0.03012243700520494,
        "ciLow": 0.3164279150082958,
        "ciHigh": 3.3597625611821815,
        "hedgesG": 1.9534108135663144,
        "pValueAdj": 0.03012243700520494
      },
      {
        "meanDiff": 3.9214285714285726,
        "statistic": 9.682374038919848,
        "df": 6.384644421604948,
        "pValueRaw": 4.753739379366951e-05,
        "ciLow": 2.944706384752858,
        "ciHigh": 4.898150758104287,
        "hedgesG": 5.319716898821222,
        "pValueAdj": 9.507478758733902e-05
      }
    ]
  },
  {
    "name": "k=5 dose-response",
    "control": [
      5.0,
      5.2,
      4.8,
      5.1,
      4.9,
      5.05
    ],
    "treated": [
      [
        9.0,
        9.5,
        8.5,
        9.2,
        8.8,
        9.1
      ],
      [
        6.0,
        6.3,
        5.7,
        6.1,
        5.9,
        6.05
      ],
      [
        5.6,
        5.9,
        5.3,
        5.7,
        5.5,
        5.65
      ],
      [
        5.3,
        5.5,
        5.1,
        5.4,
        5.2,
        5.35
      ]
    ],
    "dunnett": [
      {
        "statistic": 31.732761314901744,
        "pValueAdj": 0.0
      },
      {
        "statistic": 7.916697209538889,
        "pValueAdj": 1.630923762707681e-08
      },
      {
        "statistic": 4.750018325723334,
        "pValueAdj": 0.0002751087581422951
      },
      {
        "statistic": 2.375009162861667,
        "pValueAdj": 0.0814542544013307
      }
    ],
    "holmWelch": [
      {
        "meanDiff": 4.008333333333334,
        "statistic": 26.422190798003182,
        "df": 6.6844153359019485,
        "pValueRaw": 5.1699641365361224e-08,
        "ciLow": 3.6461501031186008,
        "ciHigh": 4.370516563548066,
        "hedgesG": 13.465246659104707,
        "pValueAdj": 2.067985654614449e-07
      },
      {
        "meanDiff": 0.9999999999999991,
        "statistic": 9.931270663228414,
        "df": 9.024555461473334,
        "pValueRaw": 3.713099413024925e-06,
        "ciLow": 0.772313200419234,
        "ciHigh": 1.227686799580764,
        "hedgesG": 5.241384207731871,
        "pValueAdj": 1.1139298239074776e-05
      },
      {
        "meanDiff": 0.5999999999999996,
        "statistic": 5.958762397937042,
        "df": 9.024555461473323,
        "pValueRaw": 0.00021073062446721745,
        "ciLow": 0.3723132004192343,
        "ciHigh": 0.827686799580765,
        "hedgesG": 3.1448305246391195,
        "pValueAdj": 0.0004214612489344349
      },
      {
        "meanDiff": 0.2999999999999998,
        "statistic": 3.6365491603879554,
        "df": 10.0,
        "pValueRaw": 0.004562791391099047,
        "ciLow": 0.11618792263911767,
        "ciHigh": 0.48381207736088194,
        "hedgesG": 1.9373622018207954,
        "pValueAdj": 0.004562791391099047
      }
    ]
  },
  {
    "name": "extreme variance ratio",
    "control": [
      50.0,
      50.1,
      49.9,
      50.05,
      49.95,
      50.02
    ],
    "treated": [
      [
        40.0,
        60.0,
        30.0,
        70.0,
        45.0
      ]
    ],
    "dunnett": [
      {
        "statistic": -0.15564146189803318,
        "pValueAdj": 0.8797985763559337
      }
    ],
    "holmWelch": [
      {
        "meanDiff": -1.0033333333333303,
        "statistic": -0.14049360530773747,
        "df": 4.000132462092936,
        "pValueRaw": 0.8950606498407925,
        "ciLow": -20.831022544995108,
        "ciHigh": 18.824355878328447,
        "hedgesG": -0.07089731950654399,
        "pValueAdj": 0.8950606498407925
      }
    ]
  },
  {
    "name": "near-equal means",
    "control": [
      10.0,
      10.5,
      9.5,
      10.2,
      9.8,
      10.1
    ],
    "treated": [
      [
        10.1,
        9.9,
        10.3,
        9.7,
        10.2,
        10.0
      ]
    ],
    "dunnett": [
      {
        "statistic": 0.1007074368138394,
        "pValueAdj": 0.9218033146613884
      }
    ],
    "holmWelch": [
      {
        "meanDiff": 0.01666666666666572,
        "statistic": 0.10070743681383938,
        "df": 8.426970577425829,
        "pValueRaw": 0.9221377426673641,
        "ciLow": -0.36162784189347863,
        "ciHigh": 0.39496117522681007,
        "hedgesG": 0.05278316845310324,
        "pValueAdj": 0.9221377426673641
      }
    ]
  }
];

export const DEMOG_FIXTURES: ControlComparisonsDemogFixture[] = [
  {
    "name": "WEIGHT by RACE vs \"Asian\" (alpha=0.05)",
    "groupCol": "RACE",
    "featureCol": "WEIGHT",
    "control": "Asian",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "Black",
        "statistic": 5.9909481364236346,
        "pValueAdj": 2.3589195041751054e-09
      },
      {
        "groupLabel": "Caucasian",
        "statistic": 4.84823356046593,
        "pValueAdj": 1.8522277730470549e-06
      },
      {
        "groupLabel": "Other",
        "statistic": 1.8428512871828744,
        "pValueAdj": 0.11941744597979687
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Black",
        "meanDiff": 16.88572186836518,
        "statistic": 6.974570277728997,
        "df": 192.71194742806887,
        "pValueRaw": 4.780345926604797e-11,
        "ciLow": 12.110580623431813,
        "ciHigh": 21.660863113298547,
        "hedgesG": 0.9210467563849912,
        "pValueAdj": 1.4341037779814391e-10
      },
      {
        "groupLabel": "Caucasian",
        "meanDiff": 11.391700215217128,
        "statistic": 6.517663328715227,
        "df": 74.63280236259,
        "pValueRaw": 7.462970845258482e-09,
        "ciLow": 7.90958406091033,
        "ciHigh": 14.873816369523926,
        "hedgesG": 0.6451756556865765,
        "pValueAdj": 1.4925941690516964e-08
      },
      {
        "groupLabel": "Other",
        "meanDiff": 4.7179096045197895,
        "statistic": 2.379713919031,
        "df": 121.07071183003481,
        "pValueRaw": 0.01888791754429268,
        "ciLow": 0.7929456424057513,
        "ciHigh": 8.642873566633828,
        "hedgesG": 0.28243620438267675,
        "pValueAdj": 0.01888791754429268
      }
    ]
  },
  {
    "name": "WEIGHT by RACE vs \"Other\" (alpha=0.05)",
    "groupCol": "RACE",
    "featureCol": "WEIGHT",
    "control": "Other",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "Asian",
        "statistic": -1.8428512871828744,
        "pValueAdj": 0.16762113931968503
      },
      {
        "groupLabel": "Black",
        "statistic": 6.408126888010933,
        "pValueAdj": 4.015712207205979e-10
      },
      {
        "groupLabel": "Caucasian",
        "statistic": 6.13796026073549,
        "pValueAdj": 3.201263032437396e-09
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Asian",
        "meanDiff": -4.7179096045197895,
        "statistic": -2.379713919031,
        "df": 121.07071183003481,
        "pValueRaw": 0.01888791754429268,
        "ciLow": -8.642873566633828,
        "ciHigh": -0.7929456424057513,
        "hedgesG": -0.28243620438267675,
        "pValueAdj": 0.01888791754429268
      },
      {
        "groupLabel": "Black",
        "meanDiff": 12.167812263845391,
        "statistic": 6.215198894823417,
        "df": 263.2690582948812,
        "pValueRaw": 1.996557256730247e-09,
        "ciLow": 8.312969732906105,
        "ciHigh": 16.02265479478468,
        "hedgesG": 0.6108381356178028,
        "pValueAdj": 3.993114513460494e-09
      },
      {
        "groupLabel": "Caucasian",
        "meanDiff": 6.673790610697338,
        "statistic": 6.587973134316383,
        "df": 410.9370275899463,
        "pValueRaw": 1.367708071824397e-10,
        "ciLow": 4.682430377801243,
        "ciHigh": 8.665150843593434,
        "hedgesG": 0.34794524474954985,
        "pValueAdj": 4.1031242154731907e-10
      }
    ]
  },
  {
    "name": "WEIGHT by RACE vs \"Asian\" (alpha=0.01)",
    "groupCol": "RACE",
    "featureCol": "WEIGHT",
    "control": "Asian",
    "alpha": 0.01,
    "dunnett": [
      {
        "groupLabel": "Black",
        "statistic": 5.9909481364236346,
        "pValueAdj": 2.3589195041751054e-09
      },
      {
        "groupLabel": "Caucasian",
        "statistic": 4.84823356046593,
        "pValueAdj": 1.8522277730470549e-06
      },
      {
        "groupLabel": "Other",
        "statistic": 1.8428512871828744,
        "pValueAdj": 0.11941744597979687
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Black",
        "meanDiff": 16.88572186836518,
        "statistic": 6.974570277728997,
        "df": 192.71194742806887,
        "pValueRaw": 4.780345926604797e-11,
        "ciLow": 10.587186361567845,
        "ciHigh": 23.184257375162517,
        "hedgesG": 0.9210467563849912,
        "pValueAdj": 1.4341037779814391e-10
      },
      {
        "groupLabel": "Caucasian",
        "meanDiff": 11.391700215217128,
        "statistic": 6.517663328715227,
        "df": 74.63280236259,
        "pValueRaw": 7.462970845258482e-09,
        "ciLow": 6.771650842442006,
        "ciHigh": 16.01174958799225,
        "hedgesG": 0.6451756556865765,
        "pValueAdj": 1.4925941690516964e-08
      },
      {
        "groupLabel": "Other",
        "meanDiff": 4.7179096045197895,
        "statistic": 2.379713919031,
        "df": 121.07071183003481,
        "pValueRaw": 0.01888791754429268,
        "ciLow": -0.470527096212507,
        "ciHigh": 9.906346305252086,
        "hedgesG": 0.28243620438267675,
        "pValueAdj": 0.01888791754429268
      }
    ]
  },
  {
    "name": "HEIGHT by RACE vs \"Caucasian\" (alpha=0.05)",
    "groupCol": "RACE",
    "featureCol": "HEIGHT",
    "control": "Caucasian",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "Asian",
        "statistic": -4.15672615792949,
        "pValueAdj": 9.84339224980868e-05
      },
      {
        "groupLabel": "Black",
        "statistic": -1.2798840559740214,
        "pValueAdj": 0.4889146841210129
      },
      {
        "groupLabel": "Other",
        "statistic": -12.964752481915504,
        "pValueAdj": 0.0
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Asian",
        "meanDiff": -5.257026016718356,
        "statistic": -3.9000976426778746,
        "df": 63.54496151988452,
        "pValueRaw": 0.00023453048003019175,
        "ciLow": -7.95018046997941,
        "ciHigh": -2.563871563457303,
        "hedgesG": -0.501632877028682,
        "pValueAdj": 0.0004690609600603835
      },
      {
        "groupLabel": "Black",
        "meanDiff": -1.0836392923316396,
        "statistic": -1.3687724715550076,
        "df": 152.62613803282488,
        "pValueRaw": 0.17308160154692523,
        "ciLow": -2.647718888204898,
        "ciHigh": 0.4804403035416187,
        "hedgesG": -0.11125966678910279,
        "pValueAdj": 0.17308160154692523
      },
      {
        "groupLabel": "Other",
        "meanDiff": -7.368024140816459,
        "statistic": -14.801335572135235,
        "df": 396.7453230301635,
        "pValueRaw": 8.841795928196054e-40,
        "ciLow": -8.346668949624876,
        "ciHigh": -6.389379332008042,
        "hedgesG": -0.7842799909915051,
        "pValueAdj": 2.652538778458816e-39
      }
    ]
  },
  {
    "name": "HEIGHT by RACE vs \"Asian\" (alpha=0.1)",
    "groupCol": "RACE",
    "featureCol": "HEIGHT",
    "control": "Asian",
    "alpha": 0.1,
    "dunnett": [
      {
        "groupLabel": "Black",
        "statistic": 2.7682920576733028,
        "pValueAdj": 0.01189924782822327
      },
      {
        "groupLabel": "Caucasian",
        "statistic": 4.15672615792949,
        "pValueAdj": 6.31581677944526e-05
      },
      {
        "groupLabel": "Other",
        "statistic": -1.5400591809009054,
        "pValueAdj": 0.21196875297847884
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Black",
        "meanDiff": 4.173386724386717,
        "statistic": 2.6943521852032877,
        "df": 105.57471136064207,
        "pValueRaw": 0.00820653981545104,
        "ciLow": 1.603054436958276,
        "ciHigh": 6.7437190118151555,
        "hedgesG": 0.41487027219399236,
        "pValueAdj": 0.01641307963090208
      },
      {
        "groupLabel": "Caucasian",
        "meanDiff": 5.257026016718356,
        "statistic": 3.9000976426778746,
        "df": 63.54496151988452,
        "pValueRaw": 0.00023453048003019175,
        "ciLow": 3.0070906014113854,
        "ciHigh": 7.506961432025326,
        "hedgesG": 0.501632877028682,
        "pValueAdj": 0.0007035914400905752
      },
      {
        "groupLabel": "Other",
        "meanDiff": -2.110998124098103,
        "statistic": -1.4852038773035596,
        "df": 78.33204637589772,
        "pValueRaw": 0.1415032881612468,
        "ciLow": -4.4768958427535654,
        "ciHigh": 0.2548995945573589,
        "hedgesG": -0.21592550881622127,
        "pValueAdj": 0.1415032881612468
      }
    ]
  },
  {
    "name": "AGE by RACE vs \"Black\" (alpha=0.05)",
    "groupCol": "RACE",
    "featureCol": "AGE",
    "control": "Black",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "Asian",
        "statistic": -2.020402523364774,
        "pValueAdj": 0.09938246934142025
      },
      {
        "groupLabel": "Caucasian",
        "statistic": -1.5112986132121624,
        "pValueAdj": 0.2719806415343411
      },
      {
        "groupLabel": "Other",
        "statistic": -0.9328426614764588,
        "pValueAdj": 0.6245621577407601
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Asian",
        "meanDiff": -3.8890658174097723,
        "statistic": -2.0958113549636908,
        "df": 135.34485652322343,
        "pValueRaw": 0.03796122228263778,
        "ciLow": -7.558860901111004,
        "ciHigh": -0.21927073370854044,
        "hedgesG": -0.2977625432297566,
        "pValueAdj": 0.11388366684791333
      },
      {
        "groupLabel": "Caucasian",
        "meanDiff": -1.6553336994201544,
        "statistic": -1.5872588350105399,
        "df": 166.53274760606126,
        "pValueRaw": 0.11435091456795772,
        "ciLow": -3.7143199843786614,
        "ciHigh": 0.4036525855383526,
        "hedgesG": -0.12466767154076346,
        "pValueAdj": 0.22870182913591544
      },
      {
        "groupLabel": "Other",
        "meanDiff": -1.209687286336326,
        "statistic": -0.9735612145314676,
        "df": 306.11998589060966,
        "pValueRaw": 0.3310431140548804,
        "ciLow": -3.6546845511551695,
        "ciHigh": 1.2353099784825177,
        "hedgesG": -0.09266433618236684,
        "pValueAdj": 0.3310431140548804
      }
    ]
  },
  {
    "name": "WEIGHT by DIS_POP vs \"RA\" (alpha=0.05)",
    "groupCol": "DIS_POP",
    "featureCol": "WEIGHT",
    "control": "RA",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "AS",
        "statistic": 2.262853583410822,
        "pValueAdj": 0.10984380830647889
      },
      {
        "groupLabel": "Indigestion",
        "statistic": -5.208283793516497,
        "pValueAdj": 9.794166402343762e-07
      },
      {
        "groupLabel": "PsA",
        "statistic": 7.780588448471501,
        "pValueAdj": 4.0190073491430667e-14
      },
      {
        "groupLabel": "Psoriasis",
        "statistic": 22.48411509657861,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "UC",
        "statistic": 3.1337736873502275,
        "pValueAdj": 0.008567969661376407
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "AS",
        "meanDiff": 2.698791895455571,
        "statistic": 2.911193032205452,
        "df": 382.2291962544831,
        "pValueRaw": 0.0038114757760876854,
        "ciLow": 0.8760556736070317,
        "ciHigh": 4.52152811730411,
        "hedgesG": 0.1647945109293362,
        "pValueAdj": 0.0038114757760876854
      },
      {
        "groupLabel": "Indigestion",
        "meanDiff": -3.8291986677376286,
        "statistic": -5.435596302550124,
        "df": 1559.3458746854708,
        "pValueRaw": 6.329150941580526e-08,
        "ciLow": -5.211001318613116,
        "ciHigh": -2.4473960168621414,
        "hedgesG": -0.21124534209287577,
        "pValueAdj": 1.898745282474158e-07
      },
      {
        "groupLabel": "PsA",
        "meanDiff": 10.770889368379756,
        "statistic": 7.921548049530171,
        "df": 230.24092349204972,
        "pValueRaw": 9.962704929399013e-14,
        "ciLow": 8.091853912269606,
        "ciHigh": 13.449924824489905,
        "hedgesG": 0.5839362813712241,
        "pValueAdj": 3.985081971759605e-13
      },
      {
        "groupLabel": "Psoriasis",
        "meanDiff": 14.801570096254025,
        "statistic": 20.290894510261868,
        "df": 2018.5531256865108,
        "pValueRaw": 1.866565728017177e-83,
        "ciLow": 13.370980141369875,
        "ciHigh": 16.232160051138173,
        "hedgesG": 0.7315944268416883,
        "pValueAdj": 9.332828640085885e-83
      },
      {
        "groupLabel": "UC",
        "meanDiff": 2.4826009068412986,
        "statistic": 3.2922015151485247,
        "df": 1193.7722249572614,
        "pValueRaw": 0.0010231658869953746,
        "ciLow": 1.0031207769339414,
        "ciHigh": 3.962081036748656,
        "hedgesG": 0.1374749937179133,
        "pValueAdj": 0.0020463317739907493
      }
    ]
  },
  {
    "name": "WEIGHT by DIS_POP vs \"AS\" (alpha=0.1)",
    "groupCol": "DIS_POP",
    "featureCol": "WEIGHT",
    "control": "AS",
    "alpha": 0.1,
    "dunnett": [
      {
        "groupLabel": "Indigestion",
        "statistic": -5.030309877818092,
        "pValueAdj": 1.2468935033371764e-06
      },
      {
        "groupLabel": "PsA",
        "statistic": 4.6145613434886075,
        "pValueAdj": 1.685300786524735e-05
      },
      {
        "groupLabel": "Psoriasis",
        "statistic": 9.637780652193674,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "RA",
        "statistic": -2.262853583410822,
        "pValueAdj": 0.07310924527856699
      },
      {
        "groupLabel": "UC",
        "statistic": -0.16244579135116727,
        "pValueAdj": 0.9997213825929049
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "Indigestion",
        "meanDiff": -6.5279905631932,
        "statistic": -6.238653354416216,
        "df": 577.2231774346984,
        "pValueRaw": 8.541622326544177e-10,
        "ciLow": -8.251896182906846,
        "ciHigh": -4.8040849434795545,
        "hedgesG": -0.40196561887189475,
        "pValueAdj": 3.416648930617671e-09
      },
      {
        "groupLabel": "PsA",
        "meanDiff": 8.072097472924185,
        "statistic": 5.1598083912964805,
        "df": 357.2995064051267,
        "pValueRaw": 4.1084792075263e-07,
        "ciLow": 5.492169555752486,
        "ciHigh": 10.652025390095883,
        "hedgesG": 0.4876911256850209,
        "pValueAdj": 1.23254376225789e-06
      },
      {
        "groupLabel": "Psoriasis",
        "meanDiff": 12.102778200798454,
        "statistic": 11.381520419059743,
        "df": 620.923681233826,
        "pValueRaw": 2.1639729875375004e-27,
        "ciLow": 10.35107516180542,
        "ciHigh": 13.854481239791488,
        "hedgesG": 0.6519194867201042,
        "pValueAdj": 1.0819864937687503e-26
      },
      {
        "groupLabel": "RA",
        "meanDiff": -2.698791895455571,
        "statistic": -2.911193032205452,
        "df": 382.2291962544831,
        "pValueRaw": 0.0038114757760876854,
        "ciLow": -4.227341425584041,
        "ciHigh": -1.170242365327102,
        "hedgesG": -0.1647945109293362,
        "pValueAdj": 0.007622951552175371
      },
      {
        "groupLabel": "UC",
        "meanDiff": -0.21619098861427233,
        "statistic": -0.20010153699637775,
        "df": 622.2191101198468,
        "pValueRaw": 0.8414665525090422,
        "ciLow": -1.9959512391085066,
        "ciHigh": 1.5635692618799606,
        "hedgesG": -0.013378681625876297,
        "pValueAdj": 0.8414665525090422
      }
    ]
  },
  {
    "name": "HEIGHT by DIS_POP vs \"Psoriasis\" (alpha=0.05)",
    "groupCol": "DIS_POP",
    "featureCol": "HEIGHT",
    "control": "Psoriasis",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "AS",
        "statistic": 3.303267327067842,
        "pValueAdj": 0.003755991277087589
      },
      {
        "groupLabel": "Indigestion",
        "statistic": -3.9909865038352343,
        "pValueAdj": 0.00026273923233577356
      },
      {
        "groupLabel": "PsA",
        "statistic": -1.5613770040953885,
        "pValueAdj": 0.3647965030796909
      },
      {
        "groupLabel": "RA",
        "statistic": -19.498647951155917,
        "pValueAdj": 0.0
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "AS",
        "meanDiff": 2.120761795558593,
        "statistic": 3.377212554452145,
        "df": 428.5598968072317,
        "pValueRaw": 0.000799093906264038,
        "ciLow": 0.886492794165376,
        "ciHigh": 3.3550307969518105,
        "hedgesG": 0.22072544260266896,
        "pValueAdj": 0.001598187812528076
      },
      {
        "groupLabel": "Indigestion",
        "meanDiff": -1.7011165160842268,
        "statistic": -3.8969406632060037,
        "df": 1902.6937734113637,
        "pValueRaw": 0.00010078770511229013,
        "ciLow": -2.5572366538407927,
        "ciHigh": -0.8449963783276606,
        "hedgesG": -0.17232596052479585,
        "pValueAdj": 0.00030236311533687036
      },
      {
        "groupLabel": "PsA",
        "meanDiff": -1.1512276643764778,
        "statistic": -1.5382450597990012,
        "df": 268.4890679049245,
        "pValueRaw": 0.12516636197890302,
        "ciLow": -2.6247131798605707,
        "ciHigh": 0.32225785110761507,
        "hedgesG": -0.11702396225106926,
        "pValueAdj": 0.12516636197890302
      },
      {
        "groupLabel": "RA",
        "meanDiff": -6.571373643480882,
        "statistic": -19.320710821287573,
        "df": 2297.8822520637227,
        "pValueRaw": 3.4388045882852243e-77,
        "ciLow": -7.23834927450643,
        "ciHigh": -5.904398012455334,
        "hedgesG": -0.680014499036136,
        "pValueAdj": 1.3755218353140897e-76
      }
    ]
  },
  {
    "name": "AGE by DIS_POP vs \"Indigestion\" (alpha=0.1)",
    "groupCol": "DIS_POP",
    "featureCol": "AGE",
    "control": "Indigestion",
    "alpha": 0.1,
    "dunnett": [
      {
        "groupLabel": "AS",
        "statistic": 2.2210919590565252,
        "pValueAdj": 0.10811502154171704
      },
      {
        "groupLabel": "PsA",
        "statistic": 9.1412201054086,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "Psoriasis",
        "statistic": 10.444084557872374,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "RA",
        "statistic": 27.546463189812027,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "UC",
        "statistic": 4.866794674209937,
        "pValueAdj": 8.38923205181974e-06
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "AS",
        "meanDiff": 1.9033592363703065,
        "statistic": 2.5924707925092285,
        "df": 523.6486703820905,
        "pValueRaw": 0.009795169683387667,
        "ciLow": 0.6935883017687858,
        "ciHigh": 3.1131301709718264,
        "hedgesG": 0.17164463075303818,
        "pValueAdj": 0.009795169683387667
      },
      {
        "groupLabel": "PsA",
        "meanDiff": 8.908305084745763,
        "statistic": 9.474756080936986,
        "df": 291.74240564507926,
        "pValueRaw": 9.61801816345373e-19,
        "ciLow": 7.3568630387988,
        "ciHigh": 10.459747130692724,
        "hedgesG": 0.7445162197428538,
        "pValueAdj": 2.8854054490361185e-18
      },
      {
        "groupLabel": "Psoriasis",
        "meanDiff": 5.750910543802831,
        "statistic": 10.73391512344559,
        "df": 1962.2804426512193,
        "pValueRaw": 3.6961921913041806e-26,
        "ciLow": 4.86923094857414,
        "ciHigh": 6.632590139031522,
        "hedgesG": 0.47259633149563085,
        "pValueAdj": 1.4784768765216722e-25
      },
      {
        "groupLabel": "RA",
        "meanDiff": 13.376736457294783,
        "statistic": 28.597341408902103,
        "df": 1616.6125960323814,
        "pValueRaw": 6.689974243952415e-146,
        "ciLow": 12.606896039784681,
        "ciHigh": 14.146576874804884,
        "hedgesG": 1.1003833220531196,
        "pValueAdj": 3.3449871219762076e-145
      },
      {
        "groupLabel": "UC",
        "meanDiff": 3.030997392438067,
        "statistic": 4.684038425531929,
        "df": 1439.103900364602,
        "pValueRaw": 3.078346989044423e-06,
        "ciLow": 1.9659424126711889,
        "ciHigh": 4.096052372204945,
        "hedgesG": 0.23598290878340436,
        "pValueAdj": 6.156693978088846e-06
      }
    ]
  },
  {
    "name": "AGE by DIS_POP vs \"PsA\" (alpha=0.05)",
    "groupCol": "DIS_POP",
    "featureCol": "AGE",
    "control": "PsA",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "AS",
        "statistic": -6.065097617959177,
        "pValueAdj": 1.9525990824220685e-09
      },
      {
        "groupLabel": "Indigestion",
        "statistic": -9.1412201054086,
        "pValueAdj": 0.0
      },
      {
        "groupLabel": "Psoriasis",
        "statistic": -3.3230649177785643,
        "pValueAdj": 0.003172793948421737
      },
      {
        "groupLabel": "RA",
        "statistic": 4.888899993577224,
        "pValueAdj": 2.349232180343108e-06
      },
      {
        "groupLabel": "UC",
        "statistic": -5.91456007313621,
        "pValueAdj": 4.625524741008746e-09
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "AS",
        "meanDiff": -7.004945848375456,
        "statistic": -6.656213728994121,
        "df": 386.0735807543602,
        "pValueRaw": 9.663673601006186e-11,
        "ciLow": -9.07408271056089,
        "ciHigh": -4.935808986190022,
        "hedgesG": -0.624245348606132,
        "pValueAdj": 3.865469440402474e-10
      },
      {
        "groupLabel": "Indigestion",
        "meanDiff": -8.908305084745763,
        "statistic": -9.474756080936986,
        "df": 291.74240564507926,
        "pValueRaw": 9.61801816345373e-19,
        "ciLow": -10.758768591451423,
        "ciHigh": -7.057841578040103,
        "hedgesG": -0.7445162197428538,
        "pValueAdj": 4.809009081726865e-18
      },
      {
        "groupLabel": "Psoriasis",
        "meanDiff": -3.1573945409429314,
        "statistic": -3.413554352415208,
        "df": 274.73012420034223,
        "pValueRaw": 0.000737874223177203,
        "ciLow": -4.97830062282723,
        "ciHigh": -1.3364884590586328,
        "hedgesG": -0.25642564980702237,
        "pValueAdj": 0.000737874223177203
      },
      {
        "groupLabel": "RA",
        "meanDiff": 4.46843137254902,
        "statistic": 5.035995544536193,
        "df": 233.73265842457252,
        "pValueRaw": 9.502378464965013e-07,
        "ciLow": 2.7203065753181903,
        "ciHigh": 6.2165561697798495,
        "hedgesG": 0.3631159728068414,
        "pValueAdj": 1.9004756929930025e-06
      },
      {
        "groupLabel": "UC",
        "meanDiff": -5.877307692307696,
        "statistic": -5.915215298566101,
        "df": 355.16155272631937,
        "pValueRaw": 7.797062700255348e-09,
        "ciLow": -7.831370163446225,
        "ciHigh": -3.9232452211691657,
        "hedgesG": -0.4529992642165935,
        "pValueAdj": 2.3391188100766045e-08
      }
    ]
  },
  {
    "name": "WEIGHT by SEX vs \"F\" (alpha=0.05)",
    "groupCol": "SEX",
    "featureCol": "WEIGHT",
    "control": "F",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "M",
        "statistic": 26.805796833154215,
        "pValueAdj": 0.0
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "M",
        "meanDiff": 13.259021288901593,
        "statistic": 26.82187600341146,
        "df": 5591.686082925005,
        "pValueRaw": 3.4946041464924556e-149,
        "ciLow": 12.289930591207652,
        "ciHigh": 14.228111986595534,
        "hedgesG": 0.7053112191400421,
        "pValueAdj": 3.4946041464924556e-149
      }
    ]
  },
  {
    "name": "HEIGHT by SEX vs \"M\" (alpha=0.05)",
    "groupCol": "SEX",
    "featureCol": "HEIGHT",
    "control": "M",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "F",
        "statistic": -64.69497647792778,
        "pValueAdj": 0.0
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "F",
        "meanDiff": -13.778390372206104,
        "statistic": -63.91083204694077,
        "df": 4443.39275878722,
        "pValueRaw": 0.0,
        "ciLow": -14.20104966340724,
        "ciHigh": -13.355731081004969,
        "hedgesG": -1.8213347896923175,
        "pValueAdj": 0.0
      }
    ]
  },
  {
    "name": "AGE by SEX vs \"F\" (alpha=0.05)",
    "groupCol": "SEX",
    "featureCol": "AGE",
    "control": "F",
    "alpha": 0.05,
    "dunnett": [
      {
        "groupLabel": "M",
        "statistic": -6.390201325424978,
        "pValueAdj": 1.786327752384409e-10
      }
    ],
    "holmWelch": [
      {
        "groupLabel": "M",
        "meanDiff": -2.266142034712466,
        "statistic": -6.412683428277025,
        "df": 5649.458572471943,
        "pValueRaw": 1.5458859894922843e-10,
        "ciLow": -2.958911076027998,
        "ciHigh": -1.573372993396934,
        "hedgesG": -0.16838307284041082,
        "pValueAdj": 1.5458859894922843e-10
      }
    ]
  }
];
