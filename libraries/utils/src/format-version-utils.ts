import * as grok from 'datagrok-api/grok';
import * as semver from 'semver';

export const DEFAULT_FLOAT_FORMAT = '#0.###';
export const G7_FLOAT_FORMAT = 'G7';

// The platform's G7 number format is supported by js-api starting from 1.27.7.
const G7_MIN_API_VERSION = '1.27.7';
let defaultFloatFormat: string | undefined;

/** Resolves the default float format used across compute inputs, outputs and exports.
 * Gated on the runtime js-api version: `G7` on platforms >= 1.27.7, `#0.###` otherwise.
 * Cached per session (the platform version does not change at runtime). */
export function getDefaultFloatFormat(): string {
  if (defaultFloatFormat === undefined) {
    const apiVersion = semver.coerce(grok.shell.build.client.version);
    defaultFloatFormat = apiVersion && semver.gte(apiVersion, G7_MIN_API_VERSION) ?
      G7_FLOAT_FORMAT : DEFAULT_FLOAT_FORMAT;
  }
  return defaultFloatFormat;
}
