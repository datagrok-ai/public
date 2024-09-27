import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {MonomerLibSummaryType} from '@datagrok-libraries/bio/src/types/index';

export const LIB_PATH = 'System:AppData/Bio/monomer-libraries/';
export const SETS_PATH: string = 'System:AppData/Bio/monomer-sets/';

export const HELM_JSON_SCHEMA_PATH = 'System:AppData/Bio/tests/libraries/HELMmonomerSchema.json';

export const LIB_SETTINGS_FOR_TESTS: UserLibSettings =
  {explicit: ['HELMCoreLibrary.json', 'polytool-lib.json'], exclude: [], duplicateMonomerPreferences: {}};

/** Summary for settings {@link LIB_SETTINGS_FOR_TESTS} */
export const monomerLibForTestsSummary: MonomerLibSummaryType = {'PEPTIDE': 326, 'RNA': 383, 'CHEM': 0};

