import * as DG from 'datagrok-api/dg';
import {HitDesignTemplate} from '../types';
import {HitDesignMolColName, PeptiHitHelmColName, TileCategoriesColName, ViDColName} from '../consts';

/** The initial one-row dataframe a new Hit Design / PeptiHit campaign starts from. */
export function createNewHitDesignDataFrame(template: HitDesignTemplate, peptiHit = false): DG.DataFrame {
  const df = DG.DataFrame.create(1);
  if (peptiHit) {
    const pepCol = df.columns.addNew(PeptiHitHelmColName, DG.TYPE.STRING);
    pepCol.semType = DG.SEMTYPE.MACROMOLECULE;
    pepCol.meta.units = 'helm';
    pepCol.setTag('units', 'helm');
    pepCol.setTag('cell.renderer', 'helm');
    pepCol.setTag('.alphabetIsMultichar', 'true');
  }

  const molCol = df.columns.addNew(HitDesignMolColName, DG.TYPE.STRING);
  molCol.semType = DG.SEMTYPE.MOLECULE;
  molCol.meta.units = DG.UNITS.Molecule.MOLBLOCK;
  molCol.set(0, null);
  if (template?.stages?.length > 0) {
    const tileCategoryCol = df.columns.addNew(TileCategoriesColName, DG.TYPE.STRING);
    tileCategoryCol.set(0, template.stages[0]);
  }
  df.columns.addNew(ViDColName, DG.TYPE.STRING);
  return df;
}
