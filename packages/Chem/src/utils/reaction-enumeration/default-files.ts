import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {_package} from '../../package';
import {detectChemSemTypes} from './shared';

export const BUNDLED_TEMPLATES = 'enumerations/reactions.csv';
export const BUNDLED_BBS = 'enumerations/bb.csv';
export const BUNDLED_EXCLUSION = 'enumerations/ex_smarts.csv';

// Chem package settings (package.json "properties", category "Enumeration"). Admins can point
// these at company-wide default files per user group via the standard package-settings mechanism.
export const TEMPLATES_PATH_SETTING = 'TemplatesPath';
export const BBS_PATH_SETTING = 'BuildingBlocksPath';
export const REAGENTS_PATH_SETTING = 'ReagentsPath';

async function csvToTable(text: string, name: string): Promise<DG.DataFrame> {
  const df = DG.DataFrame.fromCsv(text);
  df.name = name.replace(/\.csv$/i, '');
  await detectChemSemTypes(df);
  return df;
}

export async function loadBundledCsv(name: string): Promise<DG.DataFrame | null> {
  try {
    return await csvToTable(await _package.files.readAsText(name), name);
  } catch (e) {
    console.warn(`Could not load bundled file ${name}: ${e}`);
    return null;
  }
}

/** Loads the table at the admin-configured `path`; on a missing/unreadable path shows a warning
 * balloon and falls back to the bundled file (`null` for inputs with no bundled default). */
export async function loadDefaultTable(
  path: string | null | undefined, bundled: string | null, what: string,
): Promise<DG.DataFrame | null> {
  const p = path?.trim();
  if (p) {
    const fallbackNote = bundled ? ' Using the file bundled with Chem instead.' : '';
    try {
      if (await grok.dapi.files.exists(p))
        return await csvToTable(await grok.dapi.files.readAsText(p), p.split('/').pop()!);
      grok.shell.warning(`Default ${what} file not found: "${p}".${fallbackNote}`);
    } catch (e) {
      grok.shell.warning(
        `Could not load default ${what} file "${p}": ${e instanceof Error ? e.message : String(e)}.${fallbackNote}`);
    }
  }
  return bundled ? loadBundledCsv(bundled) : null;
}

export interface EnumerationDefaultTables {
  templates: DG.DataFrame | null;
  bbs: DG.DataFrame | null;
  reagents: DG.DataFrame | null;
  exclusion: DG.DataFrame | null;
}

export async function loadEnumerationDefaults(): Promise<EnumerationDefaultTables> {
  let settings: any = {};
  try {
    settings = await _package.getProperties() ?? {};
  } catch (e) {
    _package.logger.debug(`Could not read Chem settings, using bundled enumeration files: ${e}`);
  }
  const [templates, bbs, reagents, exclusion] = await Promise.all([
    loadDefaultTable(settings[TEMPLATES_PATH_SETTING], BUNDLED_TEMPLATES, 'reaction templates'),
    loadDefaultTable(settings[BBS_PATH_SETTING], BUNDLED_BBS, 'building blocks'),
    loadDefaultTable(settings[REAGENTS_PATH_SETTING], null, 'reagents'),
    loadBundledCsv(BUNDLED_EXCLUSION),
  ]);
  return {templates, bbs, reagents, exclusion};
}
