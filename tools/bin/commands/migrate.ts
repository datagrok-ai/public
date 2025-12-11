import { Project, SyntaxKind, ObjectLiteralExpression, QuoteKind, IndentationText } from 'ts-morph';
import * as path from 'path';
import { FUNC_TYPES } from 'datagrok-api/src/const';

export const toCamelCase = (str: string) =>
  str.replace(/[-_ ]+(\w)/g, (_, c) => c.toUpperCase()).replace(/^[A-Z]/, (c) => c.toLowerCase());

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
    if (!call) return;

    const expr = call.getExpression().getText();
    if (!expr.includes('grok.decorators.')) return;

    const arg = call.getArguments()[0];
    if (!arg || !arg.asKind(SyntaxKind.ObjectLiteralExpression)) return;

    const obj = arg as ObjectLiteralExpression;
    const tagsProp = obj.getProperty('tags');

    if (tagsProp) {
      const tagsArray = tagsProp.getFirstDescendantByKind(SyntaxKind.ArrayLiteralExpression);
      if (!tagsArray) return;

      const allTags = tagsArray.getElements().map((v) => v.getText().replace(/['"`]/g, ''));
      const validTags = allTags.filter((t) => Object.values(FUNC_TYPES).includes(t)).map(toCamelCase);
      const remainingTags = allTags.filter((t) => !Object.values(FUNC_TYPES).includes(t));

      if (remainingTags.length === 0)
        tagsProp.remove();
      else {
        const tagsAssignment = tagsProp.asKindOrThrow(SyntaxKind.PropertyAssignment);
        tagsAssignment.setInitializer(`[${remainingTags.map((t) => `'${t}'`).join(', ')}]`);
      }

      if (validTags.length > 0) {
        const metaProp = obj.getProperty('meta');

        if (!metaProp) {
          obj.addPropertyAssignment({
            name: 'meta',
            initializer: (writer) => {
              writer.write('{roles: [');
              validTags.forEach((t, i) => {
                if (i > 0) writer.write(', ');
                writer.write(`'${t}'`);
              });
              writer.write(']}');
            },
          });
        } else {
          const metaObj = metaProp.getFirstChildByKind(SyntaxKind.ObjectLiteralExpression);
          if (!metaObj) return;

          const otherProps = metaObj.getProperties().filter((p) => {
            return p.asKind(SyntaxKind.PropertyAssignment)?.getName() !== 'roles';
          });

          const otherText = otherProps.map((p) => p.getText()).join(', ');
          metaProp.asKindOrThrow(SyntaxKind.PropertyAssignment).setInitializer((writer) => {
            writer.write('{');
            if (otherText)
              writer.write(otherText + ', ');

            writer.write('roles: [');
            validTags.forEach((t, i) => {
              if (i > 0) writer.write(', ');
              writer.write(`'${t}'`);
            });
            writer.write(']}');
          });
        }
      }
    }
  });

  source.saveSync();
}
