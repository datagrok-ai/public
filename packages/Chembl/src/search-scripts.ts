import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {matchAndParseQuery, powerSearchQueryTable}
  from '@datagrok-libraries/db-explorer/src/search/search-widget-utils';

export async function chemblBioactivityForTargetsSearch(s: string) {
  const matchResults = matchAndParseQuery('Bioactrivity for bacterial targets for ${target}', s);
  if (!matchResults || !matchResults['target'])
    return null;

  const targetName = matchResults['target'];
  const query = DG.Func.find({package: 'Chembl', name: 'BioactivityDataForBacterialTargetsForOrganism'})[0];
  if (!query)
    return null;

  const queryCall = query.prepare({organism: targetName});
  return powerSearchQueryTable(queryCall, {
    postProcess: (tv: DG.TableView) => {
      tv.barChart({
        valueColumnName: 'standard_value', valueAggrType: 'avg', title: 'Average Standard Value',
        splitColumnName: 'standard_type', axisType: 'logarithmic', colorColumnName: 'standard_value',
      });
      tv.pieChart({categoryColumnName: 'standard_type'});
    },
  });
}

export async function chemblPKForDrugSearch(s: string) {
  const matchResults = matchAndParseQuery('Pharmacokinetic Data for ${drug}', s);
  if (!matchResults || !matchResults['drug'])
    return null;

  const drugName = matchResults['drug'];
  const query = DG.Func.find({package: 'Chembl', name: 'PKDataFromCuratedDrugPharmacokineticDataSourceForDrug'})[0];
  if (!query)
    return null;

  const queryCall = query.prepare({drug: drugName});
  return powerSearchQueryTable(queryCall, {
    postProcess: (tv: DG.TableView) => {
      tv.barChart({
        valueColumnName: 'standard_value', valueAggrType: 'count', title: 'Genders',
        splitColumnName: 'gender', colorColumnName: 'standard_value',
      });
      tv.pieChart({categoryColumnName: 'dose'});
    },

  });
}

export async function activityDetailsForTarget(s: string) {
  const matchResults = matchAndParseQuery('Compound activity details for target ${target}', s);
  if (!matchResults || !matchResults['target'])
    return null;

  const targetID = matchResults['target'];
  const query = DG.Func.find({package: 'Chembl', name: 'CompoundActivityDetailsForTarget'})[0];
  if (!query)
    return null;

  const queryCall = query.prepare({target: targetID});
  return powerSearchQueryTable(queryCall, {
    postProcess: (tv: DG.TableView) => {
    //   tv.barChart({
    //     valueColumnName: 'standard_value', valueAggrType: 'count', title: 'Genders',
    //     splitColumnName: 'gender', colorColumnName: 'standard_value',
    //   });
    //   tv.pieChart({categoryColumnName: 'dose'});
    },

  });
}
