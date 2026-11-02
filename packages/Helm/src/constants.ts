/** Names of package properties declared in the `properties` section of `package.json`. */
export const enum HelmPackagePropertiesNames {
  MonomersPerRow = 'MonomersPerRow',
}

/**
 * Fallback for `MonomersPerRow` when the package property is missing or
 * unparseable (a package loaded before its settings, an old install). Kept in
 * step with the `defaultValue` in package.json and with hwe's own default.
 */
export const DEFAULT_MONOMERS_PER_ROW = 20;
