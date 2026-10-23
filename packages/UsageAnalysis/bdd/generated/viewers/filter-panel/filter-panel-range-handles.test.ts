/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-range-handles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
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
import {filterPasses, filterPassesAll, openEmptyFilterPanel, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {dragAreaToArea, hoverArea, noErrors, readingHigher, readingHigherThanRemembered, readingLowerThanRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The range handles of a numeric filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-range-handles.feature", import.meta.url);
  test("The min and max handles narrow the range, and a reset brings every row back", {tag: ["@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(17, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
    await session.step(18, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(19, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
    await session.step(22, "When user remembers the \"min of AGE\" reading of filter panel", () => rememberReading(page, "min of AGE", el("filter panel")));
    await session.step(23, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(24, "And user drags the \"range min handle\" area of histogram viewer in filter panel to the \"bin 5\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in filter panel"), "bin 5"));
    await session.step(25, "Then the \"min of AGE\" reading of filter panel should be higher than remembered", () => readingHigherThanRemembered(page, "min of AGE", el("filter panel")));
    await session.step(26, "And the \"range min\" reading of histogram viewer in filter panel should be higher than before", () => readingHigher(page, "range min", el("histogram viewer in filter panel")));
    await session.step(27, "And 776 rows should pass the filter", () => filterPasses(page, 776));
    await session.step(28, "When user remembers the \"max of AGE\" reading of filter panel", () => rememberReading(page, "max of AGE", el("filter panel")));
    await session.step(29, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
    await session.step(30, "And user drags the \"range max handle\" area of histogram viewer in filter panel to the \"bin 15\" area", () => dragAreaToArea(page, "range max handle", el("histogram viewer in filter panel"), "bin 15"));
    await session.step(31, "Then the \"max of AGE\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "max of AGE", el("filter panel")));
    await session.step(32, "And 730 rows should pass the filter", () => filterPasses(page, 730));
    await session.step(33, "And the \"filtering of AGE\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of AGE", el("filter panel"), "true"));
    await session.step(34, "When user resets the filter", () => resetFilter(page));
    await session.step(35, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
});
