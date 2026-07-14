import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {parseProjectVocabulary} from '../publishing/publish-settings';

/**
 * Pure-parser tests for the admin-maintained `projectVocabulary` package
 * setting. The Share for Review dialog turns this list into its controlled
 * Project dropdown, so parsing (comma / newline / array, trim, dedup, empty →
 * gate) is the contract the dialog depends on.
 */
category('ProjectVocabulary', () => {
  test('comma-separated string parses to trimmed list', async () => {
    expect(JSON.stringify(parseProjectVocabulary('DMD-muscle-2026, SMA-cardiac-panel , ALS-cohort-2025')),
      JSON.stringify(['DMD-muscle-2026', 'SMA-cardiac-panel', 'ALS-cohort-2025']));
  });

  test('newline-separated string parses (admin pasted a list)', async () => {
    expect(JSON.stringify(parseProjectVocabulary('Project A\nProject B\r\nProject C')),
      JSON.stringify(['Project A', 'Project B', 'Project C']));
  });

  test('array value parses (future string_list property)', async () => {
    expect(JSON.stringify(parseProjectVocabulary(['P1', ' P2 ', 'P3'])),
      JSON.stringify(['P1', 'P2', 'P3']));
  });

  test('duplicates removed, order preserved', async () => {
    expect(JSON.stringify(parseProjectVocabulary('B, A, B, C, A')),
      JSON.stringify(['B', 'A', 'C']));
  });

  test('empty / blank / nullish → [] (dialog treats as "ask an admin")', async () => {
    expect(parseProjectVocabulary('').length, 0);
    expect(parseProjectVocabulary('   ,  , \n ').length, 0);
    expect(parseProjectVocabulary(undefined).length, 0);
    expect(parseProjectVocabulary(null).length, 0);
    expect(parseProjectVocabulary([]).length, 0);
  });
});
