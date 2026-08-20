import * as grok from 'datagrok-api/grok';
import {after, category, expect, test} from '@datagrok-libraries/test/src/test';
import {T_ALIGNMENT} from '../domain/constants';
import {
  discover, getLiveVersions, updateCuration,
} from '../service/publication-service';
import {cleanupTestPrograms, makeSavedRun, makeTestProgram, publishRun} from './fixtures';

category('ArtifactAlignment: curation', () => {
  const alignment = () => grok.dapi.domains.table(T_ALIGNMENT);

  test('publish stores tags, compounds, and path; facets and filters see them', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `t${Date.now()}`;
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Curated',
      tags: [`${marker}-pk`, `${marker}-final`], compounds: [`GRK-${marker}`], path: 'PK/exposure',
    });
    const row: any = (await getLiveVersions(result.publicationId))[0];
    expect(row.path, 'PK/exposure');
    expect(row.tags.map((t: any) => t.name).sort().join(','), `${marker}-final,${marker}-pk`);
    expect(row.compounds.map((c: any) => c.name).join(','), `GRK-${marker}`);

    const byTag = await alignment().query({filter: [
      {property: 'tags.name', operator: '=', value: `${marker}-pk`},
      {property: 'status', operator: '=', value: 'approved'},
    ]});
    expect(byTag.some((r: any) => r.id === result.rowId), true, 'tag filter must find the row');

    const byCompound = await alignment().query(
      {filter: [{property: 'compounds.registration_code', operator: '=', value: `GRK-${marker}`}]});
    expect(byCompound.some((r: any) => r.id === result.rowId), true, 'compound filter must find the row');

    const facets: any = await alignment().facets({
      filter: [{property: 'program_id', operator: '=', value: program.id}],
      facets: [
        {id: 'tags', kind: 'categories', column: 'tags.name'},
        {id: 'path', kind: 'categories', column: 'path'},
      ],
    });
    expect(facets.facets.tags.categories.some((c: any) => c.value === `${marker}-pk` || c.display === `${marker}-pk`),
      true, 'tags facet must count per tag name');
    expect(facets.facets.path.categories.some((c: any) => c.value === 'PK/exposure'), true);
  });

  test('post-publish curation edits tags and path in place', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `u${Date.now()}`;
    const result = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Retagged',
      tags: [`${marker}-old`], path: 'PK/old',
    });
    await updateCuration(result.rowId, {tags: [`${marker}-new1`, `${marker}-new2`], path: 'PK/new'});
    const row: any = (await getLiveVersions(result.publicationId))[0];
    expect(row.path, 'PK/new');
    expect(row.tags.map((t: any) => t.name).sort().join(','), `${marker}-new1,${marker}-new2`);
  });

  test('curation is copied forward to the next version on republish', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `f${Date.now()}`;
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Fwd',
      tags: [`${marker}-carried`], path: 'PK/carried', workstream: 'clinical',
    });
    await updateCuration(v1.rowId, {path: 'PK/edited'});
    const v2 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Fwd'});
    const row: any = (await getLiveVersions(v2.publicationId)).find((r: any) => r.revision === 2);
    expect(row.path, 'PK/edited', 'path must carry forward from the latest live version');
    expect(row.tags.map((t: any) => t.name).join(','), `${marker}-carried`);
    expect(row.workstream, 'clinical');
  });

  test('archive snapshots curation as plain text', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `s${Date.now()}`;
    const v1 = await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Snap',
      tags: [`${marker}-a`, `${marker}-b`], compounds: [`GRK-${marker}`],
    });
    await publishRun({sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Snap'});
    const {getPublicationHistory} = await import('../service/publication-service');
    const archived = (await getPublicationHistory(v1.publicationId))[0];
    expect(archived.tags_snapshot?.split(',').sort().join(','), `${marker}-a,${marker}-b`);
    expect(archived.compounds_snapshot, `GRK-${marker}`);
  });

  test('queryDf flattens tags and compounds with the multi-value separator annotation', async () => {
    // The platform's standard multi-value widgets (chips renderer, per-value split in
    // the categorical filter) key on the .multi-value-separator column tag; the domain
    // layer stamps it on flattened M:N columns, so Table Views get them for free.
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    const marker = `m${Date.now()}`;
    await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'MultiValue',
      tags: [`${marker}-x`, `${marker}-y`], compounds: [`GRK-${marker}`],
    });
    const df = await alignment().queryDf({
      filter: [{property: 'program_id', operator: '=', value: program.id}],
      expand: ['tags', 'compounds'],
    });
    for (const name of ['tags', 'compounds']) {
      const col = df.col(name)!;
      expect(col.meta.multiValueSeparator != null && col.meta.multiValueSeparator.trim() === ',',
        true, `'${name}' must carry the multi-value separator, got '${col.meta.multiValueSeparator}'`);
    }
    const values = `${df.col('tags')!.get(0)}`.split(',').map((s) => s.trim()).sort();
    expect(values.join(','), `${marker}-x,${marker}-y`);
  });

  test('discovery slice carries curation for facet-driven reuse pickers', async () => {
    const program = await makeTestProgram();
    const run = await makeSavedRun();
    await publishRun({
      sourceMetaCallId: run.metaCallId, programId: program.id, name: 'Slice', path: 'PK/x'});
    const slice = await discover(program.id);
    expect(slice.length, 1);
    expect(slice[0].path, 'PK/x');
    expect(Array.isArray(slice[0].tags), true, 'relations must be expanded in the slice');
  });

  after(async () => {
    await cleanupTestPrograms();
  });
}, {timeout: 120000});
