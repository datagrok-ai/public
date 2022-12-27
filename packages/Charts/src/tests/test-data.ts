import * as DG from 'datagrok-api/dg';


export const energyUK = DG.DataFrame.fromCsv(
  `value,source,target
124.72899627685547,Agricultural 'waste',Bio-conversion
0.597000002861023,Bio-conversion,Liquid
26.86199951171875,Bio-conversion,Losses
280.3219909667969,Bio-conversion,Solid
81.14399719238281,Bio-conversion,Gas
121.06600189208984,Liquid,Industry
135.8350067138672,Liquid,Road transport
3.640000104904175,Liquid,Agriculture
4.413000106811523,Liquid,Rail transport
128.69000244140625,Liquid,International shipping
14.458000183105469,Liquid,Domestic aviation
206.26699829101562,Liquid,International aviation
33.21799850463867,Liquid,National navigation
46.47700119018555,Solid,Industry
0.8820000290870667,Solid,Agriculture
400.1199951171875,Solid,Thermal generation
1.4010000228881836,Gas,Losses
48.58000183105469,Gas,Industry
0.1289999932050705,Gas,Heating and cooling - commercial
2.0959999561309814,Gas,Agriculture
151.89100646972656,Gas,Thermal generation
35.0,Biofuel imports,Liquid
35.0,Biomass imports,Solid
11.605999946594238,Coal imports,Coal
75.57099914550781,Coal,Solid
63.96500015258789,Coal reserves,Coal
10.638999938964844,District heating,Industry
22.5049991607666,District heating,Heating and cooling - commercial
46.183998107910156,District heating,Heating and cooling - homes
56.691001892089844,Electricity grid,Losses
342.1650085449219,Electricity grid,Industry
40.858001708984375,Electricity grid,Heating and cooling - commercial
113.72599792480469,Electricity grid,Heating and cooling - homes
104.4530029296875,Electricity grid,Over generation / exports
27.139999389648438,Electricity grid,H2 conversion
37.797000885009766,Electricity grid,Road transport
4.4120001792907715,Electricity grid,Agriculture
7.86299991607666,Electricity grid,Rail transport
90.00800323486328,Electricity grid,Lighting & appliances - commercial
93.49400329589844,Electricity grid,Lighting & appliances - homes
6.242000102996826,H2 conversion,Losses
20.89699935913086,H2 conversion,H2
40.71900177001953,Gas imports,Ngas
122.9520034790039,Ngas,Gas
82.23300170898438,Gas reserves,Ngas
787.1290283203125,Thermal generation,Losses
79.3290023803711,Thermal generation,District heating
525.531005859375,Thermal generation,Electricity grid
7.013000011444092,Geothermal,Electricity grid
20.89699935913086,H2,Road transport
6.994999885559082,Hydro,Electricity grid
4.375,Marine algae,Bio-conversion
839.97802734375,Nuclear,Thermal generation
504.2869873046875,Oil imports,Oil
611.989990234375,Oil,Liquid
107.7030029296875,Oil reserves,Oil
77.80999755859375,Other waste,Bio-conversion
56.58700180053711,Other waste,Solid
70.6719970703125,Pumped heat,Heating and cooling - commercial
193.0260009765625,Pumped heat,Heating and cooling - homes
59.9010009765625,Solar PV,Electricity grid
19.26300048828125,Solar Thermal,Heating and cooling - homes
59.9010009765625,Solar,Solar PV
19.26300048828125,Solar,Solar Thermal
9.45199966430664,Tidal,Electricity grid
182.00999450683594,UK land based bioenergy,Bio-conversion
19.01300048828125,Wave,Electricity grid
289.3659973144531,Wind,Electricity grid`);

export const demog = DG.DataFrame.fromCsv(
  `USUBJID,AGE,SEX,RACE,DIS_POP,HEIGHT,WEIGHT,DEMOG,CONTROL,STARTED,SEVERITY
  X0273T21000100001,53,F,Caucasian,Indigestion,160.484,73.2,53 C F,false,9/11/1990,None
  X0273T21000100002,40,F,Caucasian,Indigestion,163.646,93,40 C F,true,4/3/1991,None
  X0273T21000200001,50,F,Caucasian,Indigestion,158.027,62.7,50 C F,true,5/15/1990,High
  X0273T21000200002,30,F,Other,Indigestion,158.419,64.5,30 O F,true,2/24/1991,None
  X0273T21000300001,50,F,Caucasian,Indigestion,161.695,67.6,50 C F,false,3/5/1991,None
  X0273T21000300002,43,M,Caucasian,Indigestion,188.234,91,43 C M,false,5/23/1991,None
  X0273T21000300003,26,F,Caucasian,Indigestion,174.705,74.1,26 C F,false,8/2/1990,High
  X0273T21000300004,47,M,Caucasian,Indigestion,175.228,85.1,47 C M,true,1/21/1990,None
  X0273T21000300005,30,F,Other,Indigestion,150.288,64,30 O F,false,2/18/1991,High
  X0273T21000300006,32,M,Caucasian,Indigestion,180.298,67,32 C M,true,1/14/1991,None
  X0273T21000300007,34,M,Caucasian,Indigestion,175.351,115.3,34 C M,false,8/26/1990,Medium
  X0273T21000300008,21,F,Caucasian,Indigestion,155.091,70.5,21 C F,false,4/17/1990,Low
  X0273T21000400001,58,F,Caucasian,Indigestion,174.066,90,58 C F,true,4/5/1990,None
  X0273T21000400002,45,M,Caucasian,Indigestion,183.83,67.1,45 C M,false,7/21/1991,Medium
  X0273T21000400003,49,F,Caucasian,Indigestion,160.832,90.5,49 C F,false,7/23/1991,None
  X0273T21000400004,24,M,Caucasian,Indigestion,185.324,72,24 C M,true,4/16/1991,High
  X0273T21000400005,51,F,Caucasian,Indigestion,173.336,100.9,51 C F,false,11/27/1990,None
  X0273T21000400006,24,F,Caucasian,Indigestion,170.298,82,24 C F,false,1/28/1991,None
  X0273T21000400007,27,F,Caucasian,Indigestion,170.412,69,27 C F,false,1/27/1990,None
  X0273T21000400008,46,M,Caucasian,Indigestion,190.208,105,46 C M,false,11/17/1990,None
  X0273T21000500001,41,M,Caucasian,Indigestion,191.734,130.9,41 C M,false,8/5/1990,High
  X0273T21000500002,28,F,Caucasian,Indigestion,173.776,79.5,28 C F,false,9/4/1991,None
  X0273T21000500003,30,M,Caucasian,Indigestion,183.956,60.9,30 C M,true,3/28/1991,Medium
  X0273T21000500004,49,F,Caucasian,Indigestion,161.66,50,49 C F,true,8/25/1991,Low
  X0273T21000500005,31,F,Caucasian,Indigestion,163.741,76.4,31 C F,false,6/16/1991,None
  X0273T21000500006,51,F,Caucasian,Indigestion,164.98,60,51 C F,false,5/28/1990,None
  X0273T21000500007,73,F,Caucasian,Indigestion,162.63,52.7,73 C F,true,11/24/1990,None
  X0273T21000500008,49,M,Caucasian,Indigestion,169.074,66.4,49 C M,false,6/12/1991,High
  X0273T21000500009,28,F,Caucasian,Indigestion,168.421,48,28 C F,false,4/10/1990,None
  X0273T21000500010,37,F,Caucasian,Indigestion,163.696,58,37 C F,true,1/17/1991,None
  X0273T21000500011,28,M,Caucasian,Indigestion,163.347,65,28 C M,false,10/6/1991,None
  X0273T21000500012,52,M,Caucasian,Indigestion,175.145,77.3,52 C M,false,6/11/1991,None
  X0273T21000600001,40,F,Caucasian,Indigestion,175.372,89.9,40 C F,false,1/4/1990,None
  X0273T21000600002,53,F,Caucasian,Indigestion,163.273,62.5,53 C F,true,8/30/1991,None
  X0273T21000700001,27,F,Caucasian,Indigestion,150.081,48.2,27 C F,false,4/26/1990,Low
  X0273T21000700002,37,F,Caucasian,Indigestion,152.418,37.7,37 C F,true,1/10/1991,None
  X0273T21000700003,39,M,Caucasian,Indigestion,191.969,103.9,39 C M,false,7/10/1991,Medium
  X0273T21000700004,35,F,Caucasian,Indigestion,152.323,44.8,35 C F,true,3/8/1990,None
  X0273T21000700005,43,M,Caucasian,Indigestion,183.813,86.4,43 C M,false,11/18/1991,None
  X0273T21000700006,46,M,Caucasian,Indigestion,165.142,60.2,46 C M,false,3/22/1990,None
  X0273T21000700007,45,F,Caucasian,Indigestion,170.008,75.9,45 C F,false,3/13/1991,Low
  X0273T21000800001,66,F,Caucasian,Indigestion,159.749,84.5,66 C F,false,3/1/1991,Medium
  X0273T21000800002,45,F,Caucasian,Indigestion,163.642,50.9,45 C F,true,5/20/1991,None
  X0273T21000800003,74,M,Caucasian,Indigestion,177.848,85.1,74 C M,true,7/31/1990,None
  X0273T21000800004,25,M,Caucasian,Indigestion,172.308,65.3,25 C M,false,5/6/1990,None
  X0273T21000800005,67,F,Caucasian,Indigestion,162.789,75.7,67 C F,false,5/15/1991,Low
  X0273T21000800006,25,M,Asian,Indigestion,185.787,83,25 A M,false,9/18/1991,Low
  X0273T21000900001,46,F,Caucasian,Indigestion,168.974,118.9,46 C F,true,5/8/1991,None
  X0273T21000900002,50,F,Caucasian,Indigestion,170.259,110.5,50 C F,false,7/20/1991,None
  X0273T21000900003,28,F,Caucasian,Indigestion,175.997,59.1,28 C F,false,3/1/1990,None
  X0273T21000900004,26,F,Caucasian,Indigestion,165.032,51.3,26 C F,false,9/1/1991,None
  X0273T21000900005,40,F,Other,Indigestion,168.748,78.9,40 O F,true,6/30/1990,None
  X0273T21000900006,25,M,Caucasian,Indigestion,183.964,82.2,25 C M,false,5/10/1990,Medium
  X0273T21000900007,34,M,Caucasian,Indigestion,183.02,95.3,34 C M,false,7/1/1991,None
  X0273T21000900008,43,F,Caucasian,Indigestion,158.658,54,43 C F,false,10/14/1991,None
  X0273T21000900009,39,F,Black,Indigestion,168.439,49.3,39 B F,true,1/2/1991,None
  X0273T21000900010,42,F,Caucasian,Indigestion,165.319,68.5,42 C F,false,2/6/1990,None
  X0273T21000900011,20,F,Caucasian,Indigestion,163.398,72.1,20 C F,false,1/12/1991,None
  X0273T21000900012,52,F,Caucasian,Indigestion,160.295,83.6,52 C F,true,3/24/1990,None
  X0273T21000900013,54,M,Caucasian,Indigestion,173.158,100.9,54 C M,false,8/31/1990,Low
  X0273T21000900014,43,M,Black,Indigestion,180.98,78.2,43 B M,true,6/8/1990,High
  X0273T21000900015,49,F,Caucasian,Indigestion,166.624,72.7,49 C F,false,7/7/1990,None
  X0273T21000900016,41,F,Caucasian,Indigestion,168.327,64,41 C F,false,4/10/1991,None
  X0273T21000900017,56,M,Caucasian,Indigestion,183.848,91.7,56 C M,false,9/11/1990,None
  X0273T21000900018,27,M,Caucasian,Indigestion,175.273,67.3,27 C M,false,2/11/1991,None
  X0273T21000900019,27,F,Caucasian,Indigestion,168.654,64.1,27 C F,true,6/13/1990,None
  X0273T21000900020,53,F,Black,Indigestion,163.365,54.1,53 B F,false,1/7/1991,None
  X0273T21000900021,42,M,Caucasian,Indigestion,189.69,99.4,42 C M,false,4/12/1990,High
  X0273T21000900022,47,F,Caucasian,Indigestion,161.741,64.4,47 C F,true,8/1/1991,None
  X0273T21000900023,33,M,Caucasian,Indigestion,175.648,95,33 C M,false,7/18/1990,None
  X0273T21000900024,55,F,Caucasian,Indigestion,160.622,55.9,55 C F,false,2/28/1990,None
  X0273T21000900025,39,F,Caucasian,Indigestion,168.721,87.3,39 C F,false,2/15/1991,Low
  X0273T21000900026,22,F,Black,Indigestion,160.481,59.1,22 B F,false,2/14/1991,None
  X0273T21000900027,67,M,Caucasian,Indigestion,183.12,78.6,67 C M,false,12/18/1989,None
  X0273T21000900028,40,M,Caucasian,Indigestion,178.351,71.8,40 C M,false,2/6/1990,None
  X0273T21001000001,25,M,Caucasian,Indigestion,163.237,44.3,25 C M,false,10/29/1990,None
  X0273T21001000002,55,M,Caucasian,Indigestion,172.376,70,55 C M,false,1/8/1991,None
  X0273T21001000003,20,M,Caucasian,Indigestion,165.909,54,20 C M,false,1/6/1990,Medium
  X0273T21001000004,42,F,Caucasian,Indigestion,161.488,107.3,42 C F,false,11/3/1991,None
  X0273T21001000005,40,F,Caucasian,Indigestion,161.166,60,40 C F,false,11/28/1990,Low
  X0273T21001000006,76,F,Caucasian,Indigestion,164.632,78.6,76 C F,false,10/21/1991,None
  X0273T21001000007,23,M,Caucasian,Indigestion,169.768,73,23 C M,false,4/2/1990,None
  X0273T21001100001,23,F,Caucasian,Indigestion,159.042,57.6,23 C F,false,11/27/1990,None
  X0273T21001100002,38,M,Caucasian,Indigestion,181.167,121.7,38 C M,false,3/23/1991,None
  X0273T21001100003,27,F,Caucasian,Indigestion,165.6,52.1,27 C F,false,11/6/1990,None
  X0273T21001100004,55,M,Caucasian,Indigestion,162.376,56.7,55 C M,false,4/7/1991,Low
  X0273T21001100005,30,F,Caucasian,Indigestion,160.718,57.9,30 C F,false,11/30/1990,None
  X0273T21001100006,39,M,Caucasian,Indigestion,170.132,72.2,39 C M,false,8/23/1990,Medium
  X0273T21001100007,61,F,Caucasian,Indigestion,163.098,76.7,61 C F,false,7/24/1991,None
  X0273T21001200001,71,M,Caucasian,Indigestion,180.745,63.4,71 C M,true,8/28/1990,None
  X0273T21001200002,48,F,Caucasian,Indigestion,166.322,69.3,48 C F,false,2/26/1990,None
  X0273T21001200003,43,M,Caucasian,Indigestion,191.878,72.3,43 C M,false,12/9/1989,None
  X0273T21001200004,29,F,Caucasian,Indigestion,164.418,57.6,29 C F,false,3/19/1990,None
  X0273T21001200005,38,M,Caucasian,Indigestion,165.986,54.3,38 C M,true,7/19/1991,None
  X0273T21001200006,37,M,Caucasian,Indigestion,158.51,105.2,37 C M,false,1/24/1991,Low
  X0273T21001200007,48,M,Caucasian,Indigestion,164.288,50.5,48 C M,false,2/4/1991,None
  X0273T21001200008,48,F,Caucasian,Indigestion,158.41,45.7,48 C F,true,3/6/1990,Medium
  X0273T21001300001,41,M,Caucasian,Indigestion,173.234,75,41 C M,false,3/26/1990,Medium
  X0273T21001300002,30,F,Caucasian,Indigestion,173.131,96.4,30 C F,false,9/6/1990,Critical`);
