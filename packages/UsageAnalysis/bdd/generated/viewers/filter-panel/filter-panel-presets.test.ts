/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-presets.feature
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
import {clickOn, shouldBe, shouldContainText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPanelCount, filterPasses, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {dialogCloses, openDataset, pickSharingUser, runningAccountSignedIn, sharingUserSignedIn, signInAsSelf, signInAsSharingUser} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, clickPresetIcon, forgetSavedState, keepPresets, pickPanelMenu, presetBy, presetNamesAuthor, presetOffers, scrollPresets} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, menuDoesNotList, menuLists, menuListsAbove, noBalloons, noErrors, pickFromOpenMenu, typeIntoMenuSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {openViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter presets in the Save or Apply menu", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-presets.feature", import.meta.url);
  test("Filter presets in the Save or Apply menu", {tag: ["@journey", "@serial", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And no saved filter state \"bdd shared preset\" is kept, now or when the feature ends", () => forgetSavedState(page, "bdd shared preset"));
    await session.step(18, "And 24 saved filter presets \"bdd bulk #\" on the text column \"SEX\" are kept until the feature ends", () => keepPresets(page, 24, "bdd bulk #", "SEX"));
    await session.step(19, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(20, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await run.scenario("A preset saved with a description is listed first, as one's own", async () => {
      await session.step(23, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(24, "And user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(25, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(26, "When user picks \"Save or Apply | Save...\" from the filter panel menu", () => pickPanelMenu(page, "Save or Apply | Save..."));
      await session.step(27, "And user types \"bdd shared preset\" into Name input in \"Save filter preset\" dialog", () => typeInto(page, "bdd shared preset", el("Name input in \"Save filter preset\" dialog")));
      await session.step(28, "And user types \"Asian patients only\" into Description input in \"Save filter preset\" dialog", () => typeInto(page, "Asian patients only", el("Description input in \"Save filter preset\" dialog")));
      await session.step(29, "And user clicks on OK button in \"Save filter preset\" dialog", () => clickOn(page, el("OK button in \"Save filter preset\" dialog")));
      await session.step(30, "Then the \"Save filter preset\" dialog should close", () => dialogCloses(page, "Save filter preset"));
      await session.step(31, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(32, "Then the open menu should list \"Save or Apply > bdd shared preset\" above \"bdd bulk 24\"", () => menuListsAbove(page, "Save or Apply > bdd shared preset", "bdd bulk 24"));
      await session.step(33, "And the filter preset \"bdd shared preset\" in the open menu should be by \"me\"", () => presetBy(page, "bdd shared preset", "me"));
      await session.step(34, "And the filter preset \"bdd shared preset\" in the open menu should offer \"share and delete\"", () => presetOffers(page, "bdd shared preset", "share and delete"));
      await session.step(35, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("The menu shows the first 20 presets, the next ones on scrolling to the end, and finds any by search", async () => {
      await session.step(38, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(39, "Then the open menu should list \"Save or Apply > bdd bulk 6\"", () => menuLists(page, "Save or Apply > bdd bulk 6"));
      await session.step(40, "And the open menu should not list \"Save or Apply > bdd bulk 5\"", () => menuDoesNotList(page, "Save or Apply > bdd bulk 5"));
      await session.step(41, "When user scrolls the presets under Save or Apply to the end in the open menu", () => scrollPresets(page));
      await session.step(42, "Then the open menu should list \"Save or Apply > bdd bulk 5\"", () => menuLists(page, "Save or Apply > bdd bulk 5"));
      await session.step(43, "And the open menu should list \"Save or Apply > bdd bulk 1\"", () => menuLists(page, "Save or Apply > bdd bulk 1"));
      await session.step(44, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(45, "And user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(46, "Then the open menu should not list \"Save or Apply > bdd bulk 4\"", () => menuDoesNotList(page, "Save or Apply > bdd bulk 4"));
      await session.step(47, "When user types \"bulk 4\" into the search box of \"Save or Apply\" in the open menu", () => typeIntoMenuSearch(page, "bulk 4", "Save or Apply"));
      await session.step(48, "Then the open menu should list \"Save or Apply > bdd bulk 4\"", () => menuLists(page, "Save or Apply > bdd bulk 4"));
      await session.step(49, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("An older preset applied comes first under Sort by Used by me", async () => {
      await session.step(52, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(53, "And user types \"bulk 3\" into the search box of \"Save or Apply\" in the open menu", () => typeIntoMenuSearch(page, "bulk 3", "Save or Apply"));
      await session.step(54, "And user picks \"Save or Apply > bdd bulk 3\" from the open menu", () => pickFromOpenMenu(page, "Save or Apply > bdd bulk 3"));
      await session.step(55, "Then the filter panel should have 2 filters", () => filterPanelCount(page, 2));
      await session.step(56, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(57, "Then the open menu should list \"Save or Apply > bdd shared preset\" above \"bdd bulk 24\"", () => menuListsAbove(page, "Save or Apply > bdd shared preset", "bdd bulk 24"));
      await session.step(58, "When user picks \"Save or Apply > Sort by > Used by me\" from the open menu", () => pickFromOpenMenu(page, "Save or Apply > Sort by > Used by me"));
      await session.step(59, "Then the open menu should list \"Save or Apply > bdd bulk 3\" above \"bdd shared preset\"", () => menuListsAbove(page, "Save or Apply > bdd bulk 3", "bdd shared preset"));
      await session.step(60, "When user picks \"Save or Apply > Sort by > Updated\" from the open menu", () => pickFromOpenMenu(page, "Save or Apply > Sort by > Updated"));
      await session.step(61, "Then the open menu should not list \"Save or Apply > bdd bulk 3\"", () => menuDoesNotList(page, "Save or Apply > bdd bulk 3"));
      await session.step(62, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("A preset shared to view is applied by the second account, which can neither share nor delete it", async () => {
      await session.step(65, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(66, "And user clicks the share icon of the filter preset \"bdd shared preset\" in the open menu", () => clickPresetIcon(page, "share", "bdd shared preset"));
      await session.step(67, "Then \"Share bdd shared preset\" dialog should be visible", () => shouldBe(page, el("\"Share bdd shared preset\" dialog"), "visible"));
      await session.step(68, "When user picks the sharing user in \"User, group, or email\" input in \"Share bdd shared preset\" dialog", () => pickSharingUser(page, el("\"User, group, or email\" input in \"Share bdd shared preset\" dialog")));
      await session.step(69, "And user unchecks \"Send notifications\" input in \"Share bdd shared preset\" dialog", () => uncheck(page, el("\"Send notifications\" input in \"Share bdd shared preset\" dialog")));
      await session.step(70, "And user clicks on OK button in \"Share bdd shared preset\" dialog", () => clickOn(page, el("OK button in \"Share bdd shared preset\" dialog")));
      await session.step(71, "Then the \"Share bdd shared preset\" dialog should close", () => dialogCloses(page, "Share bdd shared preset"));
      await session.step(72, "When user signs in as the sharing user", () => signInAsSharingUser(page));
      await session.step(73, "Then the sharing user should be signed in", () => sharingUserSignedIn(page));
      await session.step(74, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(75, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(76, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(77, "Then the open menu should not list \"Save or Apply > bdd bulk 24\"", () => menuDoesNotList(page, "Save or Apply > bdd bulk 24"));
      await session.step(78, "And the filter preset \"bdd shared preset\" in the open menu should name its author", () => presetNamesAuthor(page, "bdd shared preset"));
      await session.step(79, "And the filter preset \"bdd shared preset\" in the open menu should offer \"nothing\"", () => presetOffers(page, "bdd shared preset", "nothing"));
      await session.step(80, "When user picks \"Save or Apply > bdd shared preset\" from the open menu", () => pickFromOpenMenu(page, "Save or Apply > bdd shared preset"));
      await session.step(81, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(82, "When user signs in as themselves again", () => signInAsSelf(page));
      await session.step(83, "Then the running account should be signed in", () => runningAccountSignedIn(page));
    });
    await run.scenario("Deleting one's own preset asks first, then takes it out of the menu", async () => {
      await session.step(86, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(87, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(88, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(89, "And user clicks the delete icon of the filter preset \"bdd shared preset\" in the open menu", () => clickPresetIcon(page, "delete", "bdd shared preset"));
      await session.step(90, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(91, "And \"Are you sure?\" dialog should contain text \"deleted for everyone it is shared with\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "deleted for everyone it is shared with"));
      await session.step(92, "When user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(93, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(94, "When user opens the viewer menu of filter panel", () => openViewerMenu(page, el("filter panel")));
      await session.step(95, "Then the open menu should list \"Save or Apply > bdd bulk 24\"", () => menuLists(page, "Save or Apply > bdd bulk 24"));
      await session.step(96, "And the open menu should not list \"Save or Apply > bdd shared preset\"", () => menuDoesNotList(page, "Save or Apply > bdd shared preset"));
      await session.step(97, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(98, "Then no errors should have been logged", () => noErrors(page));
      await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
