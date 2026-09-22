/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-files.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, followingShouldBe, isExpanded, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Files section of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-files.feature", import.meta.url);
  test("The Files section lists its file shares", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(32, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]), [["Files---App-Data tree node inside browse tree"],["Files---Demo tree node inside browse tree"]]);
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A folder opens as a folder view of its own", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(39, "When user clicks on Files---Demo tree node inside browse tree", () => clickOn(page, el("Files---Demo tree node inside browse tree")));
    await session.step(40, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
    await session.step(42, "And gallery should contain text \"chem\"", () => shouldContainText(page, el("gallery"), "chem"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A tabular file opens as a preview with its rows", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(47, "Given Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(48, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(49, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(50, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
    await session.step(52, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
