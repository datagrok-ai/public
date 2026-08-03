import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {withAttachedView} from '../helpers';

// CardView config mirrors (hierarchyProperties, renderMode, objectType, CustomCardView subclass).
// Getters read the rendered DOM; unattached CardView getters throw on a null context.
// A `dataSource` is required: attaching the view triggers an initial refresh, which throws
// "DataSource not specified" as an unhandled async rejection when none is supplied.
category('AI: App: CardView Extras', () => {
  const withCardView = (body: (cv: DG.CardView) => void): Promise<void> =>
    // awaitRender=false: these getters are model-backed (attachment suffices, no render needed).
    withAttachedView<DG.CardView>(() => DG.CardView.create({dataSource: grok.dapi.scripts}), (cv) => {
      expect(cv instanceof DG.CardView, true);
      body(cv);
    }, false);

  test('hierarchyProperties map round-trips', async () => {
    await withCardView((cv) => {
      cv.hierarchyProperties = {author: 'author', tag: 'tags'};
      expect(cv.hierarchyProperties != null, true);
      expect(cv.hierarchyProperties['author'], 'author');
      expect(cv.hierarchyProperties['tag'], 'tags');
    });
  });

  test('renderMode round-trips through RENDER_MODE values', async () => {
    await withCardView((cv) => {
      cv.renderMode = DG.RENDER_MODE.GRID;
      expect(cv.renderMode, DG.RENDER_MODE.GRID);
      cv.renderMode = DG.RENDER_MODE.CARD;
      expect(cv.renderMode, DG.RENDER_MODE.CARD);
      cv.renderMode = DG.RENDER_MODE.BRIEF;
      expect(cv.renderMode, DG.RENDER_MODE.BRIEF);
    });
  });

  test('objectType round-trips', async () => {
    await withCardView((cv) => {
      cv.objectType = 'Script';
      expect(cv.objectType, 'Script');
      cv.objectType = 'Project';
      expect(cv.objectType, 'Project');
    });
  });

  test('CustomCardView subclass constructs a valid CardView', async () => {
    await withAttachedView(() => new DG.CustomCardView({dataSource: grok.dapi.scripts}), async (ccv) => {
      expect(ccv instanceof DG.CustomCardView, true);
      expect(ccv instanceof DG.CardView, true);
      // Exercise a model getter to prove the attached subclass is a live CardView.
      // GRID mode would drive refreshGrid() on the construction-time async init() refresh,
      // which needs Dart-side grid properties a JS CardView has no meta for — use BRIEF.
      ccv.renderMode = DG.RENDER_MODE.BRIEF;
      expect(ccv.renderMode, DG.RENDER_MODE.BRIEF);
    }, false);
  });
}, {owner: 'agolovko@datagrok.ai'});
