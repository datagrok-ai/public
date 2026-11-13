import {_package} from '../package';
import {BioPackagePropertiesNames} from '../package-types';

/** The sequence renderers' package defaults: what a column without its own setting gets */
export type BioRendererDefaults = {
  fontSize: number,
  /** The max monomer length; 'long' (no shortening) is 50, the longest the column setting takes */
  maxMonomerLength: number,
  /** Rows with connections (BILN cyclizations, HELM bonds) draw them as lines */
  showConnections: boolean,
};

/** Reads {@link BioRendererDefaults} from the package settings at each call, so a change in the settings reaches the
 * renderers (they read it when they lay their cells out anew): `_package.properties` is a snapshot taken at init. */
export function getBioRendererDefaults(): BioRendererDefaults {
  let settings: { [name: string]: any } | null = null;
  try {
    settings = _package.settings ?? null;
  } catch {
    // a package instance not registered with the platform (the test bundle's) has no settings to read
  }
  return parseBioRendererDefaults(settings);
}

/** {@link BioRendererDefaults} from package settings: values as package.json declares them (strings, or the typed
 * values), the declared defaults for missing or unreadable ones */
export function parseBioRendererDefaults(settings: { [name: string]: any } | null): BioRendererDefaults {
  const fontSize = parseInt(String(settings?.[BioPackagePropertiesNames.FontSize]));
  const maxLengthRaw = settings?.[BioPackagePropertiesNames.MaxMonomerLength];
  const maxLength = maxLengthRaw === 'long' ? 50 : parseInt(String(maxLengthRaw));
  const show = settings?.[BioPackagePropertiesNames.ShowConnections];
  return {
    fontSize: fontSize > 0 ? fontSize : 12,
    maxMonomerLength: maxLength > 0 ? maxLength : 4,
    showConnections: show === undefined || show === null || show === '' ? true : show === true || show === 'true',
  };
}
