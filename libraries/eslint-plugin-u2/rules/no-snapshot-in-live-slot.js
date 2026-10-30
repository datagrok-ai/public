'use strict';

const {U2Imports, isFunction, sourceOf, walk} = require('../lib/u2-imports');

const SIGNAL_FACTORIES = new Set(['signal', 'computed', 'lens']);
const SIGNAL_TYPES = new Set(['Signal', 'ReadonlySignal', 'Computed']);
/** What u2 takes live when no slot table is given (check-consumer derives the full one from u2's declarations). */
const DEFAULT_SLOTS = {functions: {span: {at: [0]}, label: {at: [0]}, div: {at: [0]}}, classes: {}};

const thisKey = (node) => node.type === 'MemberExpression' && !node.computed && node.object.type === 'ThisExpression' &&
  node.property.type === 'Identifier' ? `this.${node.property.name}` : null;
const keyOf = (node) => node.type === 'Identifier' ? node.name : thisKey(node);
const propName = (node) => node.type === 'Identifier' ? node.name :
  node.type === 'Literal' && typeof node.value === 'string' ? node.value : null;

module.exports = {
  meta: {
    type: 'problem',
    docs: {description: 'A slot u2 keeps live gets a signal, not the value read once'},
    messages: {
      snapshot: '`{{read}}` reads the value once inside {{slot}}, a slot u2 keeps live, so what it feeds there ' +
        'stays as it was. Hand over the signal (`{{fix}}`), or derive one: `computed(() => …)`.',
    },
    schema: [{type: 'object'}],
  },
  create(context) {
    const slots = context.options[0] ?? DEFAULT_SLOTS;
    const functions = slots.functions ?? {};
    const classes = slots.classes ?? {};
    const source = sourceOf(context);
    const keys = source.visitorKeys;

    return {
      Program(program) {
        const imports = new U2Imports(program);
        const instances = new Map();
        const signals = new Set();
        const classOf = (type) => {
          const ref = type?.typeAnnotation;
          return ref?.type === 'TSTypeReference' ? imports.resolve(ref.typeName)?.name ?? null : null;
        };
        const learn = (key, init, type) => {
          if (key === null)
            return;
          const annotated = classOf(type);
          if (annotated !== null && SIGNAL_TYPES.has(annotated))
            signals.add(key);
          else if (annotated !== null && classes[annotated])
            instances.set(key, annotated);
          const callee = init && (init.type === 'NewExpression' || init.type === 'CallExpression') ?
            imports.resolve(init.callee) : null;
          if (callee === null)
            return;
          if (init.type === 'CallExpression' && SIGNAL_FACTORIES.has(callee.name))
            signals.add(key);
          else if (init.type === 'NewExpression' && classes[callee.name])
            instances.set(key, callee.name);
          else if (init.type === 'CallExpression' && classes[functions[callee.name]?.returns])
            instances.set(key, functions[callee.name].returns);
        };
        walk(program, keys, (n) => {
          if (n.type === 'VariableDeclarator')
            learn(keyOf(n.id), n.init, n.id.typeAnnotation);
          else if (n.type === 'PropertyDefinition' && !n.computed && n.key.type === 'Identifier')
            learn(`this.${n.key.name}`, n.value, n.typeAnnotation);
          else if (n.type === 'AssignmentExpression' && n.operator === '=')
            learn(keyOf(n.left), n.right, null);
          else if (n.type === 'Identifier' && n.typeAnnotation && n.parent && isFunction(n.parent))
            learn(n.name, null, n.typeAnnotation);
          else if (n.type === 'TSParameterProperty' && n.parameter.type === 'Identifier')
            learn(`this.${n.parameter.name}`, null, n.parameter.typeAnnotation);
          return true;
        });

        const snapshotIn = (expr) => {
          let found = null;
          walk(expr, keys, (n) => {
            if (found !== null || isFunction(n) ||
              n !== expr && (n.type === 'CallExpression' || n.type === 'NewExpression') && slotsOf(n) !== undefined)
              return false;
            if (n.type === 'CallExpression' && n.arguments.length === 0 && n.callee.type === 'MemberExpression' &&
              !n.callee.computed && n.callee.property.name === 'peek')
              found = {node: n, fix: source.getText(n.callee.object)};
            else if (n.type === 'MemberExpression' && !n.computed && n.property.name === 'value') {
              const object = n.object;
              const key = keyOf(object);
              if (object.type === 'MemberExpression' && !object.computed && object.property.name === 'value')
                found = {node: n, fix: source.getText(object)};
              else if (key !== null && signals.has(key))
                found = {node: n, fix: key};
              else if (key !== null && classes[instances.get(key)]?.valueAccessor)
                found = {node: n, fix: `${key}.signal`};
            }
            return found === null;
          });
          return found;
        };

        const slotsOf = (call) => {
          const callee = imports.resolve(call.callee);
          if (callee !== null)
            return call.type === 'NewExpression' ? classes[callee.name]?.ctor : functions[callee.name];
          const c = call.callee;
          if (call.type !== 'CallExpression' || c.type !== 'MemberExpression' || c.computed)
            return undefined;
          const receiver = keyOf(c.object);
          return receiver === null ? undefined : classes[instances.get(receiver)]?.methods?.[c.property.name];
        };

        const report = (expr, slot) => {
          const found = snapshotIn(expr);
          if (found !== null) {
            const read = source.getText(found.node);
            context.report({node: found.node, messageId: 'snapshot',
              data: {read: read.length > 40 ? `${read.slice(0, 39)}…` : read, slot, fix: found.fix}});
          }
        };

        walk(program, keys, (n) => {
          if (n.type !== 'CallExpression' && n.type !== 'NewExpression')
            return true;
          const slot = slotsOf(n);
          if (slot === undefined)
            return true;
          const name = `${n.type === 'NewExpression' ? 'new ' : ''}${source.getText(n.callee)}`;
          for (const at of slot.at ?? []) {
            if (n.arguments[at])
              report(n.arguments[at], `\`${name}(…)\``);
          }
          for (const [at, names] of Object.entries(slot.options ?? {})) {
            const arg = n.arguments[Number(at)];
            for (const p of arg?.type === 'ObjectExpression' ? arg.properties : []) {
              if (p.type === 'Property' && !p.computed && names.includes(propName(p.key)))
                report(p.value, `\`${propName(p.key)}\` of \`${name}(…)\``);
            }
          }
          return true;
        });
      },
    };
  },
};
