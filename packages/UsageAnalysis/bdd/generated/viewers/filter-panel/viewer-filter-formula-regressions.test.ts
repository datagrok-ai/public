/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/viewer-filter-formula-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, public#344, core#1141]
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
import {filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Viewer Filter formula regressions", () => {
  const session = feature(test, "features/viewers/filter-panel/viewer-filter-formula-regressions.feature", import.meta.url);
  test("A valid viewer Filter formula filters the plot and an empty one restores it (public#344)", {tag: ["@viewers", "@realizes:viewers.filters", "@realizes:public#344"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(18, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(22, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} == \"Black\"'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} == \"Black\""));
    await session.step(23, "Then scatter plot viewer should show 27 rows", () => showsRows(page, el("scatter plot viewer"), 27));
    await session.step(24, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(25, "When user sets \"Filter\" property of scatter plot viewer to \"\"", () => setProperty(page, "Filter", el("scatter plot viewer"), ""));
    await session.step(26, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
  test("An invalid viewer Filter formula is ignored and leaves the page working (public#344)", {tag: ["@viewers", "@realizes:viewers.filters", "@known-failure:34", "@realizes:public#344", "@realizes:core#1141"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(18, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
    await knownFailure(async () => {
      await session.step(31, "When user sets \"Filter\" property of scatter plot viewer to '${RACE} = \"Black\"'", () => setProperty(page, "Filter", el("scatter plot viewer"), "${RACE} = \"Black\""));
      await session.step(32, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(33, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    }, 34);
  });
});
