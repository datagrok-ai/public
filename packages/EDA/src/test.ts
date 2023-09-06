import {SampleData, OneWayAnovaTable, 
  getFactorizedData, getLevelsStat, computeOneWayAnovaTable} from './stat-tools';

const factors = ['m', 'c', 'm', 'c', 'm', 'f', 'f', 'm', 'f', 'c', 'm', 'c', 'f']; // m - 5; f - 4; c - 4
const values = [180, 100, 175, 110, 185, 160, 165, 182, 163, 105, 182, 104, 158];

const res = getLevelsStat(factors, values);
console.log(res);

console.log('=====================================================================================');

const anova = computeOneWayAnovaTable(res);
console.log(anova);

console.log('*************************************************************************************');

const cats = ['e', 'm', 'd', 'e', 'm', 'd', 'e', 'm', 'd', 'e', 'm', 'd', 'e', 'm', 'd'];
const vals = [9, 4, 1, 12, 6, 3, 4, 8, 4, 8, 2, 5, 7, 10, 2];

const resa = getFactorizedData(cats, vals);
console.log(resa);

console.log('=====================================================================================');

const anovaa = computeOneWayAnovaTable(resa);
console.log(anovaa);