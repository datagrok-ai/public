import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Molecule3DUnitsHandler} from '@datagrok-libraries/bio/src/molecule-3d';
import {Molecule3DUnits} from '@datagrok-libraries/bio/src/molecule-3d/molecule-3d-units-handler';

export const IMPORT = {
  [Molecule3DUnits.pdb]: {
    molColName: 'molecule',
    units: Molecule3DUnits.pdb,
    unitsHandlerClass: Molecule3DUnitsHandler,
  },
  [Molecule3DUnits.pdbqt]: {
    molColName: 'molecule',
    units: Molecule3DUnits.pdbqt,
    unitsHandlerClass: Molecule3DUnitsHandler,
  },
};
