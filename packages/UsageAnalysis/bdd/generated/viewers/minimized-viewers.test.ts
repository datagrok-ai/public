/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/minimized-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.minimize]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-regressions.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewerWith, clickArea, noErrors, readingIs, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Minimized viewers — minimize, preview under a filter, restore", () => {
  const session = feature(test, "features/viewers/minimized-viewers.feature", import.meta.url);
  test("Two minimized viewers preview the filtered rows, and one of them is restored", {tag: ["@viewers", "@realizes:viewers.minimize"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","HEIGHT"],["Y","WEIGHT"],["Color","RACE"]]), [["X","HEIGHT"],["Y","WEIGHT"],["Color","RACE"]]);
    await session.step(23, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(25, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(26, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(27, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(28, "And the \"rows shown\" reading of scatter plot viewer should be 5098", () => readingIs(page, "rows shown", el("scatter plot viewer"), 5098));
    await session.step(31, "When user hovers over scatter plot viewer", () => hoverOver(page, el("scatter plot viewer")));
    await session.step(32, "And user clicks scatter plot minimize icon", () => clickOn(page, el("scatter plot minimize icon")));
    await session.step(33, "Then the open tableview should have 0 scatter plot viewers", () => viewerCount(page, 0, "scatter plot"));
    await session.step(34, "And minimized scatter plot icon should be visible", () => shouldBe(page, el("minimized scatter plot icon"), "visible"));
    await session.step(35, "When user hovers over histogram viewer", () => hoverOver(page, el("histogram viewer")));
    await session.step(36, "And user clicks histogram minimize icon", () => clickOn(page, el("histogram minimize icon")));
    await session.step(37, "Then the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
    await session.step(38, "And minimized histogram icon should be visible", () => shouldBe(page, el("minimized histogram icon"), "visible"));
    await session.step(39, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
    await session.step(40, "Then 72 rows should pass the filter", () => filterPasses(page, 72));
    await session.step(41, "When user hovers over minimized scatter plot icon", () => hoverOver(page, el("minimized scatter plot icon")));
    await session.step(43, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(44, "And the \"rows shown\" reading of scatter plot viewer should be 63", () => readingIs(page, "rows shown", el("scatter plot viewer"), 63));
    await session.step(45, "When user hovers over minimized histogram icon", () => hoverOver(page, el("minimized histogram icon")));
    await session.step(46, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(47, "When user clicks minimized scatter plot icon", () => clickOn(page, el("minimized scatter plot icon")));
    await session.step(48, "Then the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(49, "And minimized scatter plot icon should be absent", () => shouldBe(page, el("minimized scatter plot icon"), "absent"));
    await session.step(50, "And minimized histogram icon should be visible", () => shouldBe(page, el("minimized histogram icon"), "visible"));
    await session.step(51, "And the \"rows shown\" reading of scatter plot viewer should be 63", () => readingIs(page, "rows shown", el("scatter plot viewer"), 63));
    await session.step(52, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
    await session.step(53, "And user clicks filter panel reset icon", () => clickOn(page, el("filter panel reset icon")));
    await session.step(54, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(55, "And the \"rows shown\" reading of scatter plot viewer should be 5098", () => readingIs(page, "rows shown", el("scatter plot viewer"), 5098));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
});
