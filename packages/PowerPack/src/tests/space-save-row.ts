import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {spaceSaveRow} from '@datagrok-libraries/u2/dg';

// What the platform's Save project dialog (xamgle project_meta.dart, through jsu.JsSpaceSaveRow) calls on the row.
category('Space save row', () => {
  /** A root space for the test; `made` collects the spaces the test creates under it, deleted first. */
  async function withSpace(body: (space: DG.Project, made: string[]) => Promise<void>): Promise<void> {
    const space = await grok.dapi.spaces.createRootSpace(`PowerPackSpaceRow${DG.Utils.randomString(6)}`);
    const made: string[] = [];
    try {
      await body(space, made);
    } finally {
      for (const id of made.reverse())
        await grok.dapi.spaces.delete(await grok.dapi.projects.find(id));
      await grok.dapi.spaces.delete(space);
    }
  }

  test('provides the spaceSaveRow role', async () => {
    const funcs = DG.Func.find({package: 'PowerPack', meta: {role: DG.FUNC_TYPES.SPACE_SAVE_ROW}});
    expect(funcs.map((f) => f.name).join(), 'spaceSaveRow');
  });

  test('a new project left out of a space stays in the user\'s projects', async () => {
    const row = spaceSaveRow(null, 'Kinase panel');
    try {
      expect(JSON.stringify(await row.resolveTarget()), JSON.stringify({spaceId: null, link: false, changed: false}));
      expect(row.problem, null);
    } finally {
      row.detach();
    }
  });

  test('a project re-saved into its own space stays; another space moves it, or links it', async () => {
    await withSpace(async (space, made) => {
      const other = await grok.dapi.spaces.id(space.id).addSubspace('Other');
      made.push(other.id);
      const row = spaceSaveRow(space, 'Kinase panel');
      try {
        expect((await row.resolveTarget()).changed, false);
        row.target = other;
        expect(JSON.stringify(await row.resolveTarget()),
          JSON.stringify({spaceId: other.id, link: false, changed: true}));
        row.mode = 'link';
        expect((await row.resolveTarget()).link, true);
      } finally {
        row.detach();
      }
    });
  });

  test('a space that holds a project of that name holds OK, a link does not', async () => {
    await withSpace(async (space) => {
      const project = DG.Project.create();
      // as the Save project dialog stores it: the friendly name, and the name derived from it
      project.friendlyName = 'Space row clash';
      project.name = 'SpaceRowClash';
      await grok.dapi.projects.save(project);
      try {
        await grok.dapi.spaces.id(space.id).addEntity(project.id);
        const row = spaceSaveRow(null, 'Space row clash');
        try {
          row.target = space;
          await waitUntil(() => row.problem !== null, 'the clash to be reported');
          expect(row.problem!.includes('already has a project named "Space row clash"'), true, row.problem!);
          row.setName('Space row clash 2');
          await waitUntil(() => row.problem === null, 'another name to be accepted');
        } finally {
          row.detach();
        }
      } finally {
        await grok.dapi.projects.delete(project);
      }
    });
  });
});
