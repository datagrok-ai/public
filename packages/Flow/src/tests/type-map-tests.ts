import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {
  areTypesCompatible, dgTypeToSlotType, getSlotColor, DG_TYPE_MAP,
  isStringListType, stringListToArrayLiteral,
  domainSection, domainCategory, isDomainOperation,
  CHEMINFORMATICS_PACKAGES, BIOINFORMATICS_PACKAGES,
  ROLE_COLORS, CATEGORY_COLORS, FUNC_NAME_COLORS, DEFAULT_NODE_COLOR,
  pastelize, categoricalColor,
} from '../types/type-map';

category('Flow: type-map', () => {
  test('identical types are compatible', async () => {
    expect(areTypesCompatible('dataframe', 'dataframe'), true);
    expect(areTypesCompatible('string', 'string'), true);
    expect(areTypesCompatible('column', 'column'), true);
  });

  test('dynamic and object are wildcards', async () => {
    expect(areTypesCompatible('dynamic', 'dataframe'), true);
    expect(areTypesCompatible('dataframe', 'dynamic'), true);
    expect(areTypesCompatible('object', 'int'), true);
    expect(areTypesCompatible('int', 'object'), true);
  });

  test('numeric widening is symmetric', async () => {
    expect(areTypesCompatible('int', 'double'), true);
    expect(areTypesCompatible('double', 'int'), true);
    expect(areTypesCompatible('num', 'int'), true);
    expect(areTypesCompatible('int', 'num'), true);
  });

  test('list and string_list interconvert', async () => {
    expect(areTypesCompatible('list', 'string_list'), true);
    expect(areTypesCompatible('string_list', 'list'), true);
  });

  test('list<string> folds to the string_list slot type', async () => {
    expect(dgTypeToSlotType('list<string>'), 'string_list');
    expect(dgTypeToSlotType('string_list'), 'string_list');
  });

  test('isStringListType matches string_list / list<string> only', async () => {
    expect(isStringListType('string_list'), true);
    expect(isStringListType('list<string>'), true);
    expect(isStringListType('list'), false, 'plain list (may be non-strings) is excluded');
    expect(isStringListType('num_list'), false);
    expect(isStringListType('column_list'), false);
  });

  test('stringListToArrayLiteral trims, drops empties, JSON-quotes', async () => {
    expect(stringListToArrayLiteral(' a, b ,, c '), '["a", "b", "c"]');
    expect(stringListToArrayLiteral(''), '[]');
    expect(stringListToArrayLiteral('   '), '[]');
    expect(stringListToArrayLiteral('x'), '["x"]');
    expect(stringListToArrayLiteral('a "q", b'), '["a \\"q\\"", "b"]', 'inner quotes escaped');
  });

  test('incompatible types are rejected', async () => {
    expect(areTypesCompatible('dataframe', 'string'), false);
    expect(areTypesCompatible('column', 'dataframe'), false);
    expect(areTypesCompatible('bool', 'datetime'), false);
    expect(areTypesCompatible('string', 'int'), false);
  });

  test('dgTypeToSlotType maps known and unknown types', async () => {
    expect(dgTypeToSlotType('dataframe'), 'dataframe');
    expect(dgTypeToSlotType('blob'), 'byte_array'); // blob slot type is byte_array
    expect(dgTypeToSlotType('totally_unknown'), 'totally_unknown'); // passthrough
  });

  test('every mapped type has a non-empty color', async () => {
    for (const [dgType, def] of Object.entries(DG_TYPE_MAP)) {
      expect(typeof def.color === 'string' && def.color.length > 0, true, `color for ${dgType}`);
      expect(getSlotColor(dgType).length > 0, true, `getSlotColor for ${dgType}`);
    }
  });

  test('pastelize lightens every channel toward white, same hue family', async () => {
    // 60% white mix of the palette hues used on title bars.
    expect(pastelize('#2ca02c'), '#abd9ab', 'input green');
    expect(pastelize('#d62728'), '#efa9a9', 'output red');
    expect(pastelize('#000000'), '#999999', 'black → 60% gray');
    expect(pastelize('#ffffff'), '#ffffff', 'white is a fixed point');
    // Every channel moves up (never down) — the pastel is strictly lighter.
    for (const c of ['#EC407A', '#7E57C2', '#4DB6AC', '#FF8A65']) {
      const orig = parseInt(c.slice(1), 16);
      const pale = parseInt(pastelize(c).slice(1), 16);
      for (const shift of [16, 8, 0])
        expect(((pale >> shift) & 0xff) >= ((orig >> shift) & 0xff), true, `${c} channel ${shift}`);
    }
    // Non-#rrggbb inputs pass through untouched.
    expect(pastelize('red'), 'red');
    expect(pastelize(''), '');
  });

  test('node identity colors come from the platform categorical palette', async () => {
    // Users must see the palette they know from every categorical coloring
    // across Datagrok (DG.Color.categoricalPalette) — no bespoke hues.
    const palette = new Set(DG.Color.categoricalPalette.map((c) => DG.Color.toHtml(c).toLowerCase()));
    expect(palette.size > 0, true, 'platform palette is available');
    const all: Array<[string, string]> = [
      ...Object.entries(ROLE_COLORS).map(([k, v]): [string, string] => [`role ${k}`, v.color]),
      ...Object.entries(CATEGORY_COLORS).map(([k, v]): [string, string] => [`category ${k}`, v.color]),
      ...Object.entries(FUNC_NAME_COLORS).map(([k, v]): [string, string] => [`func ${k}`, v.color]),
      ['default', DEFAULT_NODE_COLOR],
    ];
    for (const [what, color] of all)
      expect(palette.has(color.toLowerCase()), true, `${what} color ${color} is a palette color`);
    // categoricalColor loops over the palette like DG.Color.getCategoricalColor.
    expect(categoricalColor(0), categoricalColor(DG.Color.categoricalPalette.length));
  });

  test('domainSection routes chem/bio packages to their sections', async () => {
    expect(domainSection('Chem'), 'Cheminformatics');
    expect(domainSection('Chembl'), 'Cheminformatics');
    expect(domainSection('Admetica'), 'Cheminformatics');
    expect(domainSection('Bio'), 'Bioinformatics');
    expect(domainSection('SequenceTranslator'), 'Bioinformatics');
    expect(domainSection('BiostructureViewer'), 'Bioinformatics');
    // General packages and core get no domain (they keep their task category).
    expect(domainSection('PowerPack'), null);
    expect(domainSection('Eda'), null);
    expect(domainSection(''), null);
    expect(domainSection(null), null);
    expect(domainSection(undefined), null);
  });

  test('chem/bio package sets are disjoint and have their own colors', async () => {
    for (const p of CHEMINFORMATICS_PACKAGES)
      expect(BIOINFORMATICS_PACKAGES.has(p), false, `${p} is not in both domain sets`);
    expect(CATEGORY_COLORS['Cheminformatics'].color.length > 0, true);
    expect(CATEGORY_COLORS['Bioinformatics'].color.length > 0, true);
  });

  test('isDomainOperation requires a dataframe/column input', async () => {
    expect(isDomainOperation(['dataframe', 'string']), true);
    expect(isDomainOperation(['column']), true);
    expect(isDomainOperation(['column_list', 'double']), true);
    expect(isDomainOperation(['string', 'int']), false, 'scalar-only source is not an operation');
    expect(isDomainOperation([]), false);
  });

  test('domainCategory routes only chem/bio operations, not sources', async () => {
    // Operates on a table/column → domain section.
    expect(domainCategory('Chem', ['dataframe', 'column']), 'Cheminformatics');
    expect(domainCategory('Bio', ['column']), 'Bioinformatics');
    // A chem/bio source (scalars → table, no data input) is NOT a domain op.
    expect(domainCategory('Chem', ['string', 'int']), null, 'a molecule generator is a source');
    expect(domainCategory('Chembl', ['int']), null, 'a DB fetch is a source');
    // Non-domain packages are never routed.
    expect(domainCategory('PowerPack', ['dataframe']), null);
    expect(domainCategory('', ['dataframe']), null);
  });
});
