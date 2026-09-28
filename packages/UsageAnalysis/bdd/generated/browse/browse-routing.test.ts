/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-routing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {counterMatchesFolder, openAddress, openRememberedAddress, rememberAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {everyValueMatches} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, noProjectOnServer, openDataset, openProject, saveAsProject, urlShouldContain, urlShouldNotContain, viewHoldsViewers, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, errorBalloonText, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Opening what the address names", () => {
  const session = feature(test, "features/browse/browse-routing.feature", import.meta.url);
  test("The address of an open project opens it again after everything was closed", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "Given no project named \"bddbrowseroute{time}\" is on the server", () => noProjectOnServer(page, session.text("bddbrowseroute{time}")));
    await session.step(23, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(24, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(25, "When user saves the current view as project \"bddbrowseroute{time}\"", () => saveAsProject(page, session.text("bddbrowseroute{time}")));
    await session.step(26, "And user closes all views", () => closeAllViews(page));
    await session.step(27, "And user opens the \"bddbrowseroute{time}\" project", () => openProject(page, session.text("bddbrowseroute{time}")));
    await session.step(28, "Then the page address should contain \"bddbrowseroute\"", () => urlShouldContain(page, "bddbrowseroute"));
    await session.step(29, "When user remembers the page address", () => rememberAddress(page));
    await session.step(30, "And user closes all views", () => closeAllViews(page));
    await session.step(31, "Then the page address should not contain \"bddbrowseroute\"", () => urlShouldNotContain(page, "bddbrowseroute"));
    await session.step(32, "When user opens the remembered address", () => openRememberedAddress(page));
    await session.step(33, "Then the page address should contain \"bddbrowseroute\"", () => urlShouldContain(page, "bddbrowseroute"));
    await session.step(34, "And the current view should hold at least 2 viewers", () => viewHoldsViewers(page, 2));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
    await session.step(36, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A folder address opens that folder", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(39, "When user opens the address \"/files/System.DemoFiles/chem\"", () => openAddress(page, "/files/System.DemoFiles/chem"));
    await session.step(40, "Then the \"Demo/chem\" view should be current", () => viewIsCurrent(page, "Demo/chem"));
    await session.step(41, "And the gallery counter should show as many items as the \"System:DemoFiles/chem/\" folder holds on the server", () => counterMatchesFolder(page, "System:DemoFiles/chem/"));
    await session.step(42, "And no errors should have been logged", () => noErrors(page));
    await session.step(43, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A query address runs the query with the parameter it carries", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(46, "When user opens the address \"/func/Dbtests.PostgresByStringChoices?shipCountry=%22Germany%22\"", () => openAddress(page, "/func/Dbtests.PostgresByStringChoices?shipCountry=%22Germany%22"));
    await session.step(47, "Then the \"PostgresByStringChoices\" view should be current", () => viewIsCurrent(page, "PostgresByStringChoices"));
    await session.step(48, "And the table should have 122 rows", () => rowCount(page, 122));
    await session.step(49, "And every value of \"shipcountry\" column should match \"^Germany$\"", () => everyValueMatches(page, "shipcountry", "^Germany$"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An address that names nothing says so and leaves the shell as it was", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(54, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(55, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(56, "When user opens the address \"/p/no.such_project/none\"", () => openAddress(page, "/p/no.such_project/none"));
    await session.step(57, "Then an error balloon containing \"Unable to get project asset\" should have been shown", () => errorBalloonText(page, "Unable to get project asset"));
    await session.step(58, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
