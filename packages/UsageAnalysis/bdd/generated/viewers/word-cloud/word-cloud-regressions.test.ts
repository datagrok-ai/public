/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/word-cloud/word-cloud-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.word-cloud, core#255]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {noneSelected, onlyOfSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Word cloud regressions of closed issues", () => {
  const session = feature(test, "features/viewers/word-cloud/word-cloud-regressions.feature", import.meta.url);
  test("Clicking a word with a Map open selects its rows and logs nothing (core#255)", {tag: ["@viewers", "@realizes:viewers.word-cloud", "@realizes:core#255"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(11, "And user adds a map viewer", () => addViewer(page, "map"));
    await session.step(12, "And user adds a word cloud viewer", () => addViewer(page, "word cloud"));
    await session.step(13, "Then the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
    await session.step(14, "And the \"column\" reading of word cloud viewer should be \"MagType\"", () => readingReads(page, "column", el("word cloud viewer"), "MagType"));
    await session.step(15, "And the \"word names\" reading of word cloud viewer should be \"Mw, Mb, Ms, ML, Me, Unk\"", () => readingReads(page, "word names", el("word cloud viewer"), "Mw, Mb, Ms, ML, Me, Unk"));
    await session.step(16, "And the \"rows of word \\\"Mb\\\"\" reading of word cloud viewer should be 72", () => readingIs(page, "rows of word \"Mb\"", el("word cloud viewer"), 72));
    await session.step(17, "And no rows should be selected", () => noneSelected(page));
    await session.step(18, "When user clicks on the \"word \\\"Mb\\\"\" area of word cloud viewer", () => clickArea(page, "word \"Mb\"", el("word cloud viewer")));
    await session.step(19, "Then 72 rows should be selected", () => selectedRowCount(page, 72));
    await session.step(20, "And only rows where \"MagType\" is \"Mb\" should be selected", () => onlyOfSelected(page, "MagType", "Mb"));
    await session.step(21, "And no errors should have been logged", () => noErrors(page));
  });
});
