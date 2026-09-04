/* Custom parameter types the step expressions use. `{element}` is deliberately permissive — any
   phrase — because resolution (registry → kind → error) is the compiler's job, with a proper
   diagnostic, not the regex's. The Cucumber VS Code extension needs the same list in its
   `cucumber.parameterTypes` setting (see README). */
import {defineParameterType} from '../../src/registry.js';
import {STATES} from '../../src/runtime/assertions.js';

defineParameterType({name: 'element', regexp: /.+?/, description: 'an element phrase: "<qualifier> <kind>", a registered name, "X in Y"'});
defineParameterType({name: 'dataset', regexp: /[\w./:-]+?/, description: 'a registered dataset alias or a platform path'});
defineParameterType({name: 'viewer', regexp: /[\w -]+?/, description: 'a viewer type by its friendly name'});
defineParameterType({name: 'key', regexp: /[A-Za-z0-9+]+/, description: 'a key or chord: Enter, Escape, Control+A, ArrowDown'});
defineParameterType({name: 'state', regexp: new RegExp(STATES.join('|')), description: STATES.join(' | ')});
