/* The manifest vocabulary the editor checks while the user types — a MIRROR of the Dart
   `DomainManifestRules` (core/shared/grok_shared/lib/src/domain_manifest_rules.dart) and
   `DomainSystemColumn.reserved` (domain_types.dart), kept to the constants a rename needs. The
   server stays the authority: the dry run re-checks everything with the same values. */

export class ManifestRules {
  /** `DomainManifestRules.identifier`: schema, table and column name shape. */
  static readonly IDENTIFIER = /^[a-z][a-z0-9_]*$/;
  /** `DomainManifestRules.reservedTableNames`: shadowed by sibling `/domains/{schema}/...` routes. */
  static readonly RESERVED_TABLE_NAMES: readonly string[] = ['transaction', 'batch', 'aggregate'];
  /** `DomainManifestRules.reservedSchemaNames`: shadowed by sibling `/domains/...` routes. */
  static readonly RESERVED_SCHEMA_NAMES: readonly string[] = ['grants', 'schemas', 'filters'];
  /** `DomainSystemColumn.reserved`: the generated system columns no manifest may declare. */
  static readonly RESERVED_COLUMN_NAMES: readonly string[] = ['id', 'is_deleted', 'version', 'created_on',
    'updated_on', 'author_id', 'idempotency_key', 'data', 'tx_id'];
  /** `DomainManifestRules.reservedColumnPrefix`: the physical names of user-added extension columns. */
  static readonly RESERVED_COLUMN_PREFIX = 'x_';
  /** `DomainManifestRules.maxUserSchemaNameLength`: room for the registry prefix within 63 chars. */
  static readonly MAX_SCHEMA_NAME_LENGTH = 59;

  static checkIdentifier(name: string, what: string): string | null {
    if (name === '')
      return `${what} name is required`;
    return ManifestRules.IDENTIFIER.test(name) ? null :
      `${what} name must be lowercase letters, digits and underscores, starting with a letter`;
  }

  /** The identifier a free-text name harmonizes to: lower case, runs of anything else folded to one
   * underscore, `s_` in front when it does not start with a letter, cut to the schema length. */
  static identifier(text: string): string {
    let id = text.toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
    if (id !== '' && !/^[a-z]/.test(id))
      id = `s_${id}`;
    return id.slice(0, ManifestRules.MAX_SCHEMA_NAME_LENGTH).replace(/_+$/, '');
  }

  static checkSchemaName(name: string): string | null {
    const problem = ManifestRules.checkIdentifier(name, 'Schema');
    if (problem !== null)
      return problem;
    if (ManifestRules.RESERVED_SCHEMA_NAMES.includes(name))
      return `Schema name "${name}" is reserved for a /domains REST route`;
    return name.length > ManifestRules.MAX_SCHEMA_NAME_LENGTH ?
      `Schema name is longer than ${ManifestRules.MAX_SCHEMA_NAME_LENGTH} characters` : null;
  }

  static checkTableName(name: string): string | null {
    const problem = ManifestRules.checkIdentifier(name, 'Table');
    if (problem !== null)
      return problem;
    return ManifestRules.RESERVED_TABLE_NAMES.includes(name) ? `Table name "${name}" is reserved` : null;
  }

  static checkColumnName(name: string): string | null {
    const problem = ManifestRules.checkIdentifier(name, 'Column');
    if (problem !== null)
      return problem;
    if (ManifestRules.RESERVED_COLUMN_NAMES.includes(name))
      return `Column name "${name}" collides with a generated system column`;
    return name.startsWith(ManifestRules.RESERVED_COLUMN_PREFIX) ?
      `Column name "${name}" starts with the reserved prefix "${ManifestRules.RESERVED_COLUMN_PREFIX}"` : null;
  }
}
