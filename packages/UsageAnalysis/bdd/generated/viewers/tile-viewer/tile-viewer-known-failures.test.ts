/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer, core#1217]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewerWith} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer — the known failures", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-known-failures.feature", import.meta.url);
  test("The lane of an empty value is offered as (empty) (core#1217)", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@edge-cases", "@known-failure:26", "@realizes:core#1217"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(18, "Given user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
      await session.step(23, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","arm"]]), [["Lanes Column Name","arm"]]);
      await session.step(25, "Then the \"lane of row 1\" reading of tile viewer should be \"x\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "x"));
      await session.step(26, "And tile viewer should have a \"lane (empty)\" area", () => hasArea(page, el("tile viewer"), "lane (empty)"));
    }, 26);
  });
});
