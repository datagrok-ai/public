/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/address-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {doubleClickOn, isExpanded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress, openRememberedAddress, rememberAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeAllViews, noProjectOnServer, noSpaceOnServer, openDataset, openProject, ownProjectGone, saveAsProject, spaceOnServer, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, errorOrWarningBalloonMatching, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Address edge cases: a space, a folder, a function, a script and a project that are not there", () => {
  const session = feature(test, "features/browse/address-edge-cases.feature", import.meta.url);
  test("A space address that names nothing leads to the Spaces list", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "When user opens the address \"/s/BDDEdgeNoSpace{time}\"", () => openAddress(page, session.text("/s/BDDEdgeNoSpace{time}")));
    await session.step(29, "Then the \"Spaces\" view should be current", () => viewIsCurrent(page, "Spaces"));
    await session.step(30, "And no errors should have been logged", () => noErrors(page));
    await session.step(31, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The address of a deleted space leads to the Spaces list", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "And a space named \"BDDEdgeAddrSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(36, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(37, "When user double-clicks on BDDEdgeAddrSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeAddrSpace{time} tree node inside browse tree"))));
    await session.step(38, "Then the \"BDDEdgeAddrSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(39, "And the page address should contain \"/s/BDDEdgeAddrSpace{time}\"", () => urlShouldContain(page, session.text("/s/BDDEdgeAddrSpace{time}")));
    await session.step(40, "When user remembers the page address", () => rememberAddress(page));
    await session.step(41, "And user closes all views", () => closeAllViews(page));
    await session.step(42, "Given no space named \"BDDEdgeAddrSpace{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(43, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(44, "Then the \"Spaces\" view should be current", () => viewIsCurrent(page, "Spaces"));
    await session.step(45, "And no errors should have been logged", () => noErrors(page));
    await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A function address that names nothing says so and keeps the open view", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(49, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(50, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(51, "When user opens the address \"/func/BddNoSuchPackage.BddNoSuchFunc{time}\"", () => openAddress(page, session.text("/func/BddNoSuchPackage.BddNoSuchFunc{time}")));
    await session.step(52, "Then an error balloon containing \"Unable to find function\" should have been shown", () => errorBalloonText(page, "Unable to find function"));
    await session.step(53, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
  });
  test("The address of a deleted project says so and keeps the open view", {tag: ["@browse", "@edge-cases", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(56, "Given no project named \"bddedgeaddrproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgeaddrproj{time}")));
    await session.step(57, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(58, "When user saves the current view as project \"bddedgeaddrproj{time}\"", () => saveAsProject(page, session.text("bddedgeaddrproj{time}")));
    await session.step(59, "And user closes all views", () => closeAllViews(page));
    await session.step(60, "And user opens the \"bddedgeaddrproj{time}\" project", () => openProject(page, session.text("bddedgeaddrproj{time}")));
    await session.step(61, "Then the page address should contain \"bddedgeaddrproj{time}\"", () => urlShouldContain(page, session.text("bddedgeaddrproj{time}")));
    await session.step(62, "When user remembers the page address", () => rememberAddress(page));
    await session.step(63, "And user closes all views", () => closeAllViews(page));
    await session.step(64, "Given the user's own project \"bddedgeaddrproj{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("bddedgeaddrproj{time}")));
    await session.step(65, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(66, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(67, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(68, "Then an error balloon containing \"Unable to get project asset\" should have been shown", () => errorBalloonText(page, "Unable to get project asset"));
    await session.step(69, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
  });
  test("A folder address that names nothing says the folder is not there", {tag: ["@browse", "@edge-cases", "@realizes:views.browse", "@known-failure:80", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(78, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(79, "When user opens the address \"/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}\"", () => openAddress(page, session.text("/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}")));
      await session.step(80, "Then an error or warning balloon matching \"not found|not exist|n't exist\" should have been shown", () => errorOrWarningBalloonMatching(page, "not found|not exist|n't exist"));
      await session.step(81, "And no errors should have been logged", () => noErrors(page));
    }, 80);
  });
  test("A script address that names nothing says the script is not there", {tag: ["@browse", "@edge-cases", "@realizes:views.browse", "@known-failure:92", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(90, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(91, "When user opens the address \"/script/BddNoSuchScript{time}\"", () => openAddress(page, session.text("/script/BddNoSuchScript{time}")));
      await session.step(92, "Then an error or warning balloon matching \"not found|not exist|n't exist\" should have been shown", () => errorOrWarningBalloonMatching(page, "not found|not exist|n't exist"));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
    }, 92);
  });
});
