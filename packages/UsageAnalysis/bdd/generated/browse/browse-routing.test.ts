/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-routing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse, core#1152]
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
import {counterMatchesFolder, openAddress, openRememberedAddress, rememberAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, noProjectOnServer, noSpaceOnServer, openDataset, openProject, ownProjectGone, saveAsProject, spaceOnServer, standHasReachableConnection, urlShouldContain, urlShouldNotContain, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, errorBalloonText, errorOrWarningBalloonMatching, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening what the address names", () => {
  const session = feature(test, "features/browse/browse-routing.feature", import.meta.url);
  test("The address of an open project opens it again after everything was closed", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "Given no project named \"bddbrowseroute\" is on the server", () => noProjectOnServer(page, "bddbrowseroute"));
    await session.step(35, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(36, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(37, "When user saves the current view as project \"bddbrowseroute\"", () => saveAsProject(page, "bddbrowseroute"));
    await session.step(38, "And user closes all views", () => closeAllViews(page));
    await session.step(39, "And user opens the \"bddbrowseroute\" project", () => openProject(page, "bddbrowseroute"));
    await session.step(40, "Then the page address should contain \"bddbrowseroute\"", () => urlShouldContain(page, "bddbrowseroute"));
    await session.step(41, "When user remembers the page address", () => rememberAddress(page));
    await session.step(42, "And user closes all views", () => closeAllViews(page));
    await session.step(43, "Then the page address should not contain \"bddbrowseroute\"", () => urlShouldNotContain(page, "bddbrowseroute"));
    await session.step(44, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(45, "Then the page address should contain \"bddbrowseroute\"", () => urlShouldContain(page, "bddbrowseroute"));
    await session.step(46, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
    await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The address of a deleted project says so and keeps the open view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "Given no project named \"bddedgeaddrproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgeaddrproj{time}")));
    await session.step(52, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(53, "When user saves the current view as project \"bddedgeaddrproj{time}\"", () => saveAsProject(page, session.text("bddedgeaddrproj{time}")));
    await session.step(54, "And user closes all views", () => closeAllViews(page));
    await session.step(55, "And user opens the \"bddedgeaddrproj{time}\" project", () => openProject(page, session.text("bddedgeaddrproj{time}")));
    await session.step(56, "Then the page address should contain \"bddedgeaddrproj{time}\"", () => urlShouldContain(page, session.text("bddedgeaddrproj{time}")));
    await session.step(57, "When user remembers the page address", () => rememberAddress(page));
    await session.step(58, "And user closes all views", () => closeAllViews(page));
    await session.step(59, "Given the user's own project \"bddedgeaddrproj{time}\" is removed now and at feature end", () => ownProjectGone(page, session.text("bddedgeaddrproj{time}")));
    await session.step(60, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(61, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(62, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(63, "Then an error balloon containing \"Unable to get project asset\" should have been shown", () => errorBalloonText(page, "Unable to get project asset"));
    await session.step(64, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("A folder address opens that folder, and a space address that names nothing or a deleted space leads to the Spaces list", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(68, "When user opens the address \"/files/System.DemoFiles/chem\"", () => openAddress(page, "/files/System.DemoFiles/chem"));
    await session.step(69, "Then the \"Demo/chem\" view should be current", () => viewIsCurrent(page, "Demo/chem"));
    await session.step(70, "And the gallery counter should show as many items as the \"System:DemoFiles/chem/\" folder holds on the server", () => counterMatchesFolder(page, "System:DemoFiles/chem/"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
    await session.step(72, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(73, "When user opens the address \"/s/BDDEdgeNoSpace{time}\"", () => openAddress(page, session.text("/s/BDDEdgeNoSpace{time}")));
    await session.step(74, "Then the \"Spaces\" view should be current", () => viewIsCurrent(page, "Spaces"));
    await session.step(75, "And no errors should have been logged", () => noErrors(page));
    await session.step(76, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(77, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(78, "And a space named \"BDDEdgeAddrSpace{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(79, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(80, "When user double-clicks on BDDEdgeAddrSpace{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeAddrSpace{time} tree node inside browse tree"))));
    await session.step(81, "Then the \"BDDEdgeAddrSpace{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(82, "And the page address should contain \"/s/BDDEdgeAddrSpace{time}\"", () => urlShouldContain(page, session.text("/s/BDDEdgeAddrSpace{time}")));
    await session.step(83, "When user remembers the page address", () => rememberAddress(page));
    await session.step(84, "And user closes all views", () => closeAllViews(page));
    await session.step(85, "Given no space named \"BDDEdgeAddrSpace{time}\" is on the server", () => noSpaceOnServer(page, session.text("BDDEdgeAddrSpace{time}")));
    await session.step(86, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(87, "Then the \"Spaces\" view should be current", () => viewIsCurrent(page, "Spaces"));
    await session.step(88, "And no errors should have been logged", () => noErrors(page));
    await session.step(89, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A query address runs the query with the parameter it carries", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(93, "Given the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(94, "When user opens the address \"/func/Dbtests.PostgresByStringChoices?shipCountry=%22Germany%22\"", () => openAddress(page, "/func/Dbtests.PostgresByStringChoices?shipCountry=%22Germany%22"));
    await session.step(95, "Then the \"PostgresByStringChoices\" view should be current", () => viewIsCurrent(page, "PostgresByStringChoices"));
    await session.step(96, "And the table should have 122 rows", () => rowCount(page, 122));
    await session.step(97, "And every value of \"shipcountry\" column should match \"^Germany$\"", () => everyValueMatches(page, "shipcountry", "^Germany$"));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
    await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An address that names no project or function says so and keeps the open view, and one that names no app opens the Apps view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(102, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(103, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(104, "When user opens the address \"/p/no.such_project/none\"", () => openAddress(page, "/p/no.such_project/none"));
    await session.step(105, "Then an error balloon containing \"Unable to get project asset\" should have been shown", () => errorBalloonText(page, "Unable to get project asset"));
    await session.step(106, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(107, "And no errors should have been logged", () => noErrors(page));
    await session.step(108, "When user opens the address \"/func/BddNoSuchPackage.BddNoSuchFunc{time}\"", () => openAddress(page, session.text("/func/BddNoSuchPackage.BddNoSuchFunc{time}")));
    await session.step(109, "Then an error balloon containing \"Unable to find function\" should have been shown", () => errorBalloonText(page, "Unable to find function"));
    await session.step(110, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(111, "When user opens the address \"/apps/BddNoSuchApp{time}\"", () => openAddress(page, session.text("/apps/BddNoSuchApp{time}")));
    await session.step(112, "Then an error balloon containing \"Application not found\" should have been shown", () => errorBalloonText(page, "Application not found"));
    await session.step(113, "And the \"Apps\" view should be current", () => viewIsCurrent(page, "Apps"));
  });
  test("A folder address that names nothing says the folder is not there [kind=folder, address=/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}]", {tag: ["@browse", "@realizes:views.browse", "@known-failure:123", "@realizes:core#1152"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(121, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(122, "When user opens the address \"/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}\"", () => openAddress(page, session.text("/files/System.AppData/UsageAnalysis/bdd-edge-no-folder-{time}")));
      await session.step(123, "Then an error or warning balloon matching \"not found|not exist|n't exist\" should have been shown", () => errorOrWarningBalloonMatching(page, "not found|not exist|n't exist"));
    }, 123);
  });
  test("A script address that names nothing says the script is not there [kind=script, address=/script/BddNoSuchScript{time}]", {tag: ["@browse", "@realizes:views.browse", "@known-failure:123", "@realizes:core#1152"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await knownFailure(async () => {
      await session.step(121, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(122, "When user opens the address \"/script/BddNoSuchScript{time}\"", () => openAddress(page, session.text("/script/BddNoSuchScript{time}")));
      await session.step(123, "Then an error or warning balloon matching \"not found|not exist|n't exist\" should have been shown", () => errorOrWarningBalloonMatching(page, "not found|not exist|n't exist"));
    }, 123);
  });
});
