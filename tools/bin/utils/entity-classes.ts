/// `@grok.decorators.entity` classes → manifest tables (grok schema generate). Reads the
/// decorators statically with ts-morph, the way `migrate.ts` reads `@func`; pure — takes a
/// Project, returns plain manifest objects.
import {Node, Project, PropertyDeclaration, SyntaxKind} from 'ts-morph';

export const DOMAIN_SCHEMA_ID = 'https://datagrok.ai/schemas/domain-schema.schema.json';

export interface EntityColumn {
  name: string;
  /** The `@column` options as written; `refClass` carries a `ref: () => Class` target. */
  options: {[key: string]: any};
}

export interface EntityClass {
  className: string;
  file: string;
  schema: string;
  table: string;
  /** The `@entity` options minus `schema` and `table`. */
  options: {[key: string]: any};
  columns: EntityColumn[];
}

const ourDecorator = /^(?:grok\.)?decorators\.(entity|column)$/;
/** Manifest table keys, in the order the manifest lists them. */
const tableKeys = ['securityMode', 'promotion', 'defaultRowVisibility', 'delegate', 'softDelete', 'audit',
  'extensible', 'idempotency', 'ginIndex', 'businessKey', 'friendlyName', 'description', 'singularName',
  'pluralName', 'filters', 'relations', 'schemas'];
const columnKeys = ['type', 'ref', 'onDelete', 'required', 'unique', 'isName', 'min', 'max', 'choices', 'semType',
  'friendlyName', 'description', 'default', 'format'];
const columnTypes = ['string', 'int', 'float', 'bool', 'datetime', 'string_list', 'ref', 'user', 'group', 'file'];
const inferredTypes: {[tsType: string]: string} = {
  string: 'string', number: 'float', boolean: 'bool', Date: 'datetime', 'string[]': 'string_list',
  'Array<string>': 'string_list',
};

/** `PlateWell` → `plate_well`, `HTTPReader` → `http_reader`. */
export function snakeCase(name: string): string {
  return name.replace(/([a-z0-9])([A-Z])/g, '$1_$2').replace(/([A-Z]+)([A-Z][a-z])/g, '$1_$2').toLowerCase();
}

function propertyName(p: Node): string {
  return (p as any).getName().replace(/^['"`]|['"`]$/g, '');
}

/** The value of a literal expression; anything else is an error naming [where]. */
function literal(expr: Node, where: string): any {
  if (Node.isStringLiteral(expr) || Node.isNoSubstitutionTemplateLiteral(expr) || Node.isNumericLiteral(expr))
    return expr.getLiteralValue();
  if (Node.isPrefixUnaryExpression(expr) && expr.getOperatorToken() === SyntaxKind.MinusToken) {
    const operand = expr.getOperand();
    if (Node.isNumericLiteral(operand))
      return -operand.getLiteralValue();
  }
  if (expr.getKind() === SyntaxKind.TrueKeyword)
    return true;
  if (expr.getKind() === SyntaxKind.FalseKeyword)
    return false;
  if (Node.isArrayLiteralExpression(expr))
    return expr.getElements().map((e) => literal(e, where));
  if (Node.isObjectLiteralExpression(expr)) {
    const o: {[key: string]: any} = {};
    for (const p of expr.getProperties()) {
      if (!Node.isPropertyAssignment(p))
        throw new Error(`${where} must be a literal (found '${p.getText()}')`);
      const key = propertyName(p);
      o[key] = literal(p.getInitializerOrThrow(), `${where}.${key}`);
    }
    return o;
  }
  throw new Error(`${where} must be a string, number, boolean, array or object literal (found '${expr.getText()}')`);
}

/** The decorator's options object as plain values; [where] names the class (and property). */
function decoratorOptions(args: Node[], where: string, decorator: string, required: boolean): {[key: string]: any} {
  const arg = args[0];
  if (arg == null) {
    if (required)
      throw new Error(`${where}: @${decorator} needs an options object`);
    return {};
  }
  if (!Node.isObjectLiteralExpression(arg))
    throw new Error(`${where}: @${decorator} options must be an object literal (found '${arg.getText()}')`);
  const o: {[key: string]: any} = {};
  for (const p of arg.getProperties()) {
    if (!Node.isPropertyAssignment(p))
      throw new Error(`${where}: @${decorator} option '${p.getText()}' must be a literal`);
    const key = propertyName(p);
    const init = p.getInitializerOrThrow();
    // `ref: () => SomeClass` names an @entity class; resolved once every class is known.
    if (decorator === 'column' && key === 'ref' && Node.isArrowFunction(init)) {
      const body = init.getBody();
      if (!Node.isIdentifier(body))
        throw new Error(`${where}: @column option 'ref' must be a table name or '() => EntityClass'`);
      o.refClass = body.getText();
      continue;
    }
    o[key] = literal(init, `${where}: @${decorator} option '${key}'`);
  }
  return o;
}

function inferType(prop: PropertyDeclaration): {type?: string, text: string} {
  const text = (prop.getTypeNode()?.getText() ?? prop.getType().getNonNullableType().getText()).replace(/\s+/g, '');
  return {type: inferredTypes[text], text};
}

/** Every `@entity` class of [project], in file-path then declaration order. */
export function readEntityClasses(project: Project): EntityClass[] {
  const classes: EntityClass[] = [];
  const files = [...project.getSourceFiles()].sort((a, b) => (a.getFilePath() < b.getFilePath() ? -1 : 1));
  for (const file of files) {
    for (const cls of file.getClasses()) {
      const decorator = cls.getDecorators().find((d) => ourDecorator.exec(d.getFullName())?.[1] === 'entity');
      if (decorator == null)
        continue;
      const className = cls.getName() ?? '(anonymous class)';
      const options = decoratorOptions(decorator.getArguments(), className, 'entity', true);
      if (typeof options.schema !== 'string' || options.schema === '')
        throw new Error(`${className}: @entity needs a 'schema'`);
      if (options.table != null && typeof options.table !== 'string')
        throw new Error(`${className}: @entity option 'table' must be a string`);
      const tableOptions: {[key: string]: any} = {};
      for (const k of Object.keys(options)) {
        if (k === 'schema' || k === 'table')
          continue;
        if (!tableKeys.includes(k))
          throw new Error(`${className}: unknown @entity option '${k}'`);
        tableOptions[k] = options[k];
      }
      const columns: EntityColumn[] = [];
      for (const prop of cls.getProperties()) {
        const cd = prop.getDecorators().find((d) => ourDecorator.exec(d.getFullName())?.[1] === 'column');
        if (cd == null)
          continue;
        const where = `${className}.${prop.getName()}`;
        const raw = decoratorOptions(cd.getArguments(), where, 'column', false);
        for (const k of Object.keys(raw)) {
          if (k !== 'name' && k !== 'refClass' && !columnKeys.includes(k))
            throw new Error(`${where}: unknown @column option '${k}'`);
        }
        if (raw.type != null && !columnTypes.includes(raw.type))
          throw new Error(`${where}: column type '${raw.type}' is not one of ${columnTypes.join(', ')}`);
        if (raw.type == null) {
          const inferred = inferType(prop);
          if (inferred.type == null)
            throw new Error(`${where}: cannot infer a column type from '${inferred.text}'; declare {type: ...}`);
          raw.type = inferred.type;
        }
        columns.push({name: raw.name ?? prop.getName(), options: raw});
      }
      classes.push({className, file: file.getFilePath(), schema: options.schema,
        table: options.table ?? snakeCase(className), options: tableOptions, columns});
    }
  }
  return classes;
}

function manifestColumn(cls: EntityClass, col: EntityColumn, byClass: Map<string, EntityClass>): {[key: string]: any} {
  const where = `${cls.className}.${col.name}`;
  const o = {...col.options};
  if (o.refClass != null) {
    const target = byClass.get(o.refClass);
    if (target == null)
      throw new Error(`${where}: ref target '${o.refClass}' is not an @entity class`);
    if (target.schema !== cls.schema)
      throw new Error(`${where}: ref target '${o.refClass}' belongs to schema '${target.schema}', not '${cls.schema}'`);
    o.ref = target.table;
  }
  if (o.ref != null && o.type !== 'user' && o.type !== 'group')
    o.type = 'ref';
  if (o.type === 'ref' && o.ref == null)
    throw new Error(`${where}: a 'ref' column needs its target ('ref')`);
  if (o.type !== 'ref')
    delete o.ref;
  const column: {[key: string]: any} = {};
  for (const k of columnKeys) {
    if (o[k] !== undefined)
      column[k] = o[k];
  }
  return column;
}

/** One manifest per schema from [classes]; `version` is the package version. */
export function buildManifests(classes: EntityClass[], version: string): {[schema: string]: any} {
  const byTable = new Map<string, EntityClass>();
  for (const c of classes) {
    const other = byTable.get(c.table);
    if (other != null)
      throw new Error(`${c.className} and ${other.className} both map to table '${c.table}'`);
    byTable.set(c.table, c);
  }
  const byClass = new Map(classes.map((c) => [c.className, c]));
  const manifests: {[schema: string]: any} = {};
  for (const c of classes) {
    manifests[c.schema] ??= {$schema: DOMAIN_SCHEMA_ID, name: c.schema, version, tables: {}};
    const columns: {[name: string]: any} = {};
    for (const col of c.columns) {
      if (columns[col.name] != null)
        throw new Error(`${c.className}: two properties map to column '${col.name}'`);
      columns[col.name] = manifestColumn(c, col, byClass);
    }
    const declared = (name: string, what: string) => {
      if (columns[name] == null)
        throw new Error(`${c.className}: ${what} '${name}' is not a @column of the class`);
    };
    if (c.options.delegate != null)
      declared(String(c.options.delegate), 'delegate');
    for (const k of c.options.businessKey ?? [])
      declared(String(k), 'businessKey column');
    const named = Object.keys(columns).filter((n) => columns[n].isName === true);
    if (named.length > 1)
      throw new Error(`${c.className}: more than one isName column (${named.join(', ')})`);
    if (named.length === 1 && columns[named[0]].type !== 'string')
      throw new Error(`${c.className}: isName column '${named[0]}' must be a string`);
    const table: {[key: string]: any} = {};
    for (const k of tableKeys) {
      if (c.options[k] !== undefined)
        table[k] = c.options[k];
    }
    table.columns = columns;
    manifests[c.schema].tables[c.table] = table;
  }
  return manifests;
}
