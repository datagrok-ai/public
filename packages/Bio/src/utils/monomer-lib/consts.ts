import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {MonomerLibSummaryType} from '@datagrok-libraries/bio/src/types/monomer-library';

export const LIB_SETTINGS_FOR_TESTS: UserLibSettings =
  {explicit: ['HELMCoreLibrary.json', 'polytool-lib.json'], exclude: [], duplicateMonomerPreferences: {}};

/** Summary for settings {@link LIB_SETTINGS_FOR_TESTS} */
export const monomerLibForTestsSummary: MonomerLibSummaryType = {'PEPTIDE': 334, 'RNA': 383, 'CHEM': 0};

