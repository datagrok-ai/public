import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {
  applySmartPathwayFilter, GENERIC_PARENT_TERMS, SPECIFIC_CHILD_TERMS,
  GostResult,
} from '../analysis/enrichment';

/** Build a minimally-shaped GostResult — only the fields the smart filter
 * reads are meaningful; the rest are filler so the type lines up. */
function gost(opts: Partial<GostResult> & {name: string; source: string; p_value: number}): GostResult {
  return {
    native: opts.native ?? `NS:${opts.name}`,
    name: opts.name,
    source: opts.source,
    p_value: opts.p_value,
    significant: opts.significant ?? true,
    term_size: opts.term_size ?? 100,
    query_size: opts.query_size ?? 50,
    intersection_size: opts.intersection_size ?? 5,
    effective_domain_size: opts.effective_domain_size ?? 20000,
    precision: opts.precision ?? 0.1,
    recall: opts.recall ?? 0.05,
    intersections: opts.intersections ?? [],
  };
}

category('Proteomics: 14-05', () => {
  test('smartPathwayFilter constants match CK-omics literal lists', async () => {
    expect(GENERIC_PARENT_TERMS.length, 6);
    expect(GENERIC_PARENT_TERMS.includes('localization'), true);
    expect(GENERIC_PARENT_TERMS.includes('cellular component organization'), true);
    expect(GENERIC_PARENT_TERMS.includes('transport'), true);
    expect(GENERIC_PARENT_TERMS.includes('cellular process'), true);
    expect(GENERIC_PARENT_TERMS.includes('biological process'), true);
    expect(GENERIC_PARENT_TERMS.includes('metabolic process'), true);

    expect(SPECIFIC_CHILD_TERMS.length, 4);
    expect(SPECIFIC_CHILD_TERMS.includes('actin'), true);
    expect(SPECIFIC_CHILD_TERMS.includes('vesicle'), true);
    expect(SPECIFIC_CHILD_TERMS.includes('endocytosis'), true);
    expect(SPECIFIC_CHILD_TERMS.includes('cytoskeleton'), true);
  });

  test('smartFilterDropsGenericParentWhenChildKept', async () => {
    // NB: substring matching is literal per CK-omics. A name containing
    // 'transport' is classified as generic even if it also contains 'vesicle'
    // — pick fixture names that don't accidentally collide with a generic
    // substring when the intent is "specific child kept".
    const results: GostResult[] = [
      gost({name: 'actin cytoskeleton organization', source: 'GO:BP', p_value: 0.001}),
      gost({name: 'cellular component organization', source: 'GO:BP', p_value: 0.002}),
      gost({name: 'vesicle docking', source: 'GO:BP', p_value: 0.003}),
      gost({name: 'transport', source: 'GO:BP', p_value: 0.004}),
      gost({name: 'protein phosphorylation', source: 'GO:BP', p_value: 0.005}),
    ];
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    // Specific children + unrelated kept; the 2 generic-parents dropped.
    expect(kept.length, 3);
    expect(stats.total, 5);
    expect(stats.kept, 3);
    expect(stats.droppedParents, 2);
    expect(stats.cappedAtN, 15);
    const names = kept.map((r) => r.name);
    expect(names.includes('actin cytoskeleton organization'), true);
    expect(names.includes('vesicle docking'), true);
    expect(names.includes('protein phosphorylation'), true);
    expect(names.includes('cellular component organization'), false);
    expect(names.includes('transport'), false);
  });

  test('smartFilterKeepsGenericParentWhenNoChildFirst', async () => {
    // CK-omics behavior: generic parent at the lowest p is evaluated first,
    // when the kept list is still empty → kept unconditionally. Later child
    // does not retroactively drop it.
    const results: GostResult[] = [
      gost({name: 'biological process', source: 'GO:BP', p_value: 0.001}),
      gost({name: 'actin filament organization', source: 'GO:BP', p_value: 0.002}),
      gost({name: 'protein folding', source: 'GO:BP', p_value: 0.003}),
    ];
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    expect(kept.length, 3);
    expect(stats.droppedParents, 0);
  });

  test('smartFilterCapsGoBpAt15', async () => {
    const results: GostResult[] = [];
    for (let i = 0; i < 30; i++)
      results.push(gost({name: `pathway ${i}`, source: 'GO:BP', p_value: 0.001 * (i + 1)}));
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    expect(kept.length, 15);
    expect(stats.total, 30);
    expect(stats.kept, 15);
    expect(stats.droppedParents, 0);
    // The 15 lowest-p rows are kept.
    for (const r of kept) {
      const idx = parseInt(r.name.replace('pathway ', ''));
      expect(idx < 15, true);
    }
  });

  test('smartFilterCombinedCapForNonGoBp (Assumption A4)', async () => {
    // 5 GO:BP (all kept because <15) + 20 non-GO:BP mixed across KEGG/REAC/WP.
    // CK-omics line 4731: combined .head(maxPerSource), NOT per-source.
    const results: GostResult[] = [];
    for (let i = 0; i < 5; i++)
      results.push(gost({name: `bp ${i}`, source: 'GO:BP', p_value: 0.001 * (i + 1)}));
    const otherSources = ['KEGG', 'REAC', 'WP'];
    for (let i = 0; i < 20; i++) {
      results.push(gost({
        name: `other ${i}`, source: otherSources[i % 3], p_value: 0.01 * (i + 1),
      }));
    }
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    // 5 GO:BP + 15 non-GO:BP combined cap = 20 total
    expect(kept.length, 20);
    expect(stats.total, 25);
    expect(stats.kept, 20);
    const nonBpKept = kept.filter((r) => r.source !== 'GO:BP');
    expect(nonBpKept.length, 15);
  });

  test('smartFilterPreservesPValueSortOrder', async () => {
    const results: GostResult[] = [
      gost({name: 'a', source: 'GO:BP', p_value: 0.005}),
      gost({name: 'b', source: 'KEGG', p_value: 0.001}),
      gost({name: 'c', source: 'GO:BP', p_value: 0.003}),
      gost({name: 'd', source: 'REAC', p_value: 0.002}),
    ];
    const {kept} = applySmartPathwayFilter(results, 15);
    for (let i = 1; i < kept.length; i++)
      expect(kept[i - 1].p_value <= kept[i].p_value, true);
  });

  test('smartFilterEmptyInput', async () => {
    const {kept, stats} = applySmartPathwayFilter([], 15);
    expect(kept.length, 0);
    expect(stats.total, 0);
    expect(stats.kept, 0);
    expect(stats.droppedParents, 0);
    expect(stats.cappedAtN, 15);
  });

  test('smartFilterStatsCount', async () => {
    const results: GostResult[] = [
      gost({name: 'actin organization', source: 'GO:BP', p_value: 0.001}),
      gost({name: 'biological process', source: 'GO:BP', p_value: 0.002}),
      gost({name: 'metabolic process', source: 'GO:BP', p_value: 0.003}),
    ];
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    expect(stats.total, results.length);
    expect(stats.kept, kept.length);
    expect(stats.droppedParents, 2);
    expect(stats.cappedAtN, 15);
  });

  test('smartFilter empty GO:BP, only other sources', async () => {
    const results: GostResult[] = [];
    for (let i = 0; i < 25; i++)
      results.push(gost({name: `kegg ${i}`, source: 'KEGG', p_value: 0.001 * (i + 1)}));
    const {kept, stats} = applySmartPathwayFilter(results, 15);
    // Combined cap on non-GO:BP only.
    expect(kept.length, 15);
    expect(stats.droppedParents, 0);
  });

  test('smartFilter non-GO:BP path also sorts by p ASC', async () => {
    // Out-of-order non-GO:BP inputs; cap should keep the lowest-p subset.
    const results: GostResult[] = [
      gost({name: 'kegg_high', source: 'KEGG', p_value: 0.1}),
      gost({name: 'kegg_low', source: 'KEGG', p_value: 0.001}),
      gost({name: 'reac_mid', source: 'REAC', p_value: 0.05}),
    ];
    const {kept} = applySmartPathwayFilter(results, 2);
    // Cap=2 → keep the 2 lowest-p (kegg_low p=0.001, reac_mid p=0.05).
    expect(kept.length, 2);
    expect(kept[0].name, 'kegg_low');
    expect(kept[1].name, 'reac_mid');
  });
});
