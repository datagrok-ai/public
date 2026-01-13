import {Project, SyntaxKind, ObjectLiteralExpression, QuoteKind, IndentationText} from 'ts-morph';
import * as path from 'path';
import {FUNC_TYPES} from '../utils/const';

const FUNC_TYPE_VALUES = new Set(Object.values(FUNC_TYPES));

export const toCamelCase = (str: string) =>
  str.replace(/[-_ ]+(\w)/g, (_, c) => c.toUpperCase()).replace(/^[A-Z]/, (c) => c.toLowerCase());

function getProp(obj: ObjectLiteralExpression, name: string) {
  return obj.getProperties().find((p) => {
    if (![SyntaxKind.PropertyAssignment, SyntaxKind.ShorthandPropertyAssignment].includes(p.getKind()))
      return false;

    const n = (p as any).getName?.();
    return !!n && n.replace(/['"`]/g, '') === name;
  });
}

const parseRoleString = (value: string) =>
  value.split(',').map((v) => v.trim()).filter(Boolean);

export function migrate(argv?: string[]) {
  const FILE_PATH = path.resolve(process.cwd(), 'src/package.ts');

  const project = new Project({
    manipulationSettings: {
      quoteKind: QuoteKind.Single,
      indentationText: IndentationText.TwoSpaces,
      useTrailingCommas: true,
    },
  });

  const source = project.addSourceFileAtPath(FILE_PATH);

  source.getDescendantsOfKind(SyntaxKind.Decorator).forEach((decorator) => {
    const call = decorator.getCallExpression();
    if (!call || !call.getExpression().getText().includes('grok.decorators.'))
      return;

    const arg = call.getArguments()[0]?.asKind(SyntaxKind.ObjectLiteralExpression);
    if (!arg)
      return;

    const tagsProp = getProp(arg, 'tags');
    if (!tagsProp)
      return;

    const tagsArray = tagsProp.getFirstDescendantByKind(SyntaxKind.ArrayLiteralExpression);
    if (!tagsArray)
      return;

    const validTags: string[] = [];
    const remainingTags: string[] = [];

    tagsArray.getElements().forEach((el) => {
      const tag = el.getText().replace(/['"`]/g, '');
      FUNC_TYPE_VALUES.has(tag) ? validTags.push(toCamelCase(tag)) : remainingTags.push(tag);
    });

    if (remainingTags.length === 0)
      tagsProp.remove();
    else {
      tagsProp
        .asKindOrThrow(SyntaxKind.PropertyAssignment)
        .setInitializer(`[${remainingTags.map((t) => `'${t}'`).join(', ')}]`);
    }

    if (validTags.length === 0)
      return;

    const metaProp = getProp(arg, 'meta');
    if (!metaProp) {
      arg.addPropertyAssignment({name: 'meta', initializer: `{role: '${validTags.join(',')}'}`});
      return;
    }

    const metaObj = metaProp.getFirstChildByKind(SyntaxKind.ObjectLiteralExpression);
    if (!metaObj)
      return;

    const roleProp = getProp(metaObj, 'role');
    const existingRoles = roleProp ?
      parseRoleString(
        roleProp
          .asKindOrThrow(SyntaxKind.PropertyAssignment)
          .getInitializer()
          ?.getText()
          .replace(/['"`]/g, '') ?? '',
      ) :
      [];

    const mergedRoles = Array.from(new Set([...existingRoles, ...validTags]));

    const otherProps = metaObj
      .getProperties()
      .filter((p) => p.asKind(SyntaxKind.PropertyAssignment)?.getName() !== 'role')
      .map((p) => p.getText())
      .join(', ');

    metaProp.asKindOrThrow(SyntaxKind.PropertyAssignment).setInitializer((writer) => {
      writer.write('{');
      if (otherProps) writer.write(otherProps + ', ');
      writer.write(`role: '${mergedRoles.join(',')}'`);
      writer.write('}');
    });
  });

  source.saveSync();
}
