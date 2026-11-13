import {defineParameterType} from '../../src/registry.js';
import {STATES} from '../../src/states.js';

// the VS Code extension's `cucumber.parameterTypes` setting needs the same list (README); `{element}` takes any
// phrase on purpose: resolving it is the compiler's job, with a proper diagnostic, not the regex's
defineParameterType({name: 'element', regexp: /.+?/, description: 'an element phrase: "<qualifier> <kind>", a registered name, "X in Y"'});
defineParameterType({name: 'dataset', regexp: /[\w./:-]+?/, description: 'a registered dataset alias or a platform path'});
defineParameterType({name: 'viewer', regexp: /[\w -]+?/, description: 'a viewer type by its friendly name'});
defineParameterType({name: 'widget', regexp: /(?:.+? (?:viewer|widget)|(?:the )?grid|(?:the )?filter panel)(?: in .+?)?/,
  description: 'an element phrase that names a viewer or a widget ("box plot viewer", "results widget"), the platform\'s grid ("grid") or its filter panel ("filter panel") — steps that need one accept nothing else'});
defineParameterType({name: 'key', regexp: /[A-Za-z0-9+]+/, description: 'a key or chord: Enter, Escape, Control+A, ArrowDown; Control/Ctrl means Command on macOS, Delete/Del means Backspace there; ControlLeft and ForwardDelete send physical keys'});
defineParameterType({name: 'state', regexp: new RegExp(STATES.join('|')), description: STATES.join(' | ')});
