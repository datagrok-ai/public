/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/apps/demos-data-access-compute-curves.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {makeRowCurrentInTable, tableFilterCount, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {customFiredWith, listenCustom} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, currentViewType, packageInstalled, shellSettingIs, shellSettingPutBack, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Data Access, Compute and Curves demos open from Browse > Apps > Demo with their content", () => {
  const session = feature(test, "features/apps/demos-data-access-compute-curves.feature", import.meta.url);
  test("The Table Linking demo opens with its second grid viewer [package=Tutorials, group=Data-Access, section=Data Access, demo=Table Linking, node=Table-Linking, content=second grid viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(23, "And Apps---Demo---Data-Access tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Data-Access tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Data-Access---Table-Linking tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Data-Access---Table-Linking tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Data Access | Table Linking\"", () => customFiredWith(page, "demo-loaded", "path", "Data Access | Table Linking"));
    await session.step(27, "And the \"Table Linking\" view should be current", () => viewIsCurrent(page, "Table Linking"));
    await session.step(28, "And second grid viewer should be visible", () => shouldBe(page, el("second grid viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Diff Studio demo opens with its line chart viewer [package=DiffStudio, group=Compute, section=Compute, demo=Diff Studio, node=Diff-Studio, content=line chart viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(23, "And Apps---Demo---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Compute tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Compute---Diff-Studio tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Compute---Diff-Studio tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Compute | Diff Studio\"", () => customFiredWith(page, "demo-loaded", "path", "Compute | Diff Studio"));
    await session.step(27, "And the \"Diff Studio\" view should be current", () => viewIsCurrent(page, "Diff Studio"));
    await session.step(28, "And line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The PK-PD Modeling demo opens with its line chart viewer [package=DiffStudio, group=Compute, section=Compute, demo=PK-PD Modeling, node=PK-PD-Modeling, content=line chart viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(23, "And Apps---Demo---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Compute tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Compute---PK-PD-Modeling tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Compute---PK-PD-Modeling tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Compute | PK-PD Modeling\"", () => customFiredWith(page, "demo-loaded", "path", "Compute | PK-PD Modeling"));
    await session.step(27, "And the \"PK-PD Modeling\" view should be current", () => viewIsCurrent(page, "PK-PD Modeling"));
    await session.step(28, "And line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Bioreactor demo opens with its line chart viewer [package=DiffStudio, group=Compute, section=Compute, demo=Bioreactor, node=Bioreactor, content=line chart viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"DiffStudio\" package is installed", () => packageInstalled(page, "DiffStudio"));
    await session.step(23, "And Apps---Demo---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Compute tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Compute---Bioreactor tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Compute---Bioreactor tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Compute | Bioreactor\"", () => customFiredWith(page, "demo-loaded", "path", "Compute | Bioreactor"));
    await session.step(27, "And the \"Bioreactor\" view should be current", () => viewIsCurrent(page, "Bioreactor"));
    await session.step(28, "And line chart viewer should be visible", () => shouldBe(page, el("line chart viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Multivariate Analysis demo opens with its bar chart viewer [package=Eda, group=Compute, section=Compute, demo=Multivariate Analysis, node=Multivariate-Analysis, content=bar chart viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"Eda\" package is installed", () => packageInstalled(page, "Eda"));
    await session.step(23, "And Apps---Demo---Compute tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Compute tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Compute---Multivariate-Analysis tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Compute---Multivariate-Analysis tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Compute | Multivariate Analysis\"", () => customFiredWith(page, "demo-loaded", "path", "Compute | Multivariate Analysis"));
    await session.step(27, "And the \"Multivariate Analysis\" view should be current", () => viewIsCurrent(page, "Multivariate Analysis"));
    await session.step(28, "And bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Curve Fitting demo opens with its grid [package=Curves, group=Curves, section=Curves, demo=Curve Fitting, node=Curve-Fitting, content=grid]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(23, "And Apps---Demo---Curves tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Curves tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Curves---Curve-Fitting tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Curves---Curve-Fitting tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Curves | Curve Fitting\"", () => customFiredWith(page, "demo-loaded", "path", "Curves | Curve Fitting"));
    await session.step(27, "And the \"Curve Fitting\" view should be current", () => viewIsCurrent(page, "Curve Fitting"));
    await session.step(28, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Assay Curves demo opens with its MultiCurveViewer viewer [package=Curves, group=Curves, section=Curves, demo=Assay Curves, node=Assay-Curves, content=MultiCurveViewer viewer]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(22, "Given the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(23, "And Apps---Demo---Curves tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Curves tree node inside browse tree")));
    await session.step(24, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(25, "When user clicks on Apps---Demo---Curves---Assay-Curves tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Curves---Assay-Curves tree node inside browse tree")));
    await session.step(26, "Then the \"demo-loaded\" custom event should have fired with path \"Curves | Assay Curves\"", () => customFiredWith(page, "demo-loaded", "path", "Curves | Assay Curves"));
    await session.step(27, "And the \"Assay Curves\" view should be current", () => viewIsCurrent(page, "Assay Curves"));
    await session.step(28, "And MultiCurveViewer viewer should be visible", () => shouldBe(page, el("MultiCurveViewer viewer"), "visible"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Files demo opens the platform's files view [demo=Files, type=files]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(43, "Given Apps---Demo---Data-Access tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Data-Access tree node inside browse tree")));
    await session.step(44, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(45, "When user clicks on Apps---Demo---Data-Access---Files tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Data-Access---Files tree node inside browse tree")));
    await session.step(46, "Then the \"demo-loaded\" custom event should have fired with path \"Data Access | Files\"", () => customFiredWith(page, "demo-loaded", "path", "Data Access | Files"));
    await session.step(47, "And the \"Files\" view should be current", () => viewIsCurrent(page, "Files"));
    await session.step(48, "And the current view should be a files view", () => currentViewType(page, "files"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
    await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Databases demo opens the platform's databases view [demo=Databases, type=databases]", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(43, "Given Apps---Demo---Data-Access tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Data-Access tree node inside browse tree")));
    await session.step(44, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(45, "When user clicks on Apps---Demo---Data-Access---Databases tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Data-Access---Databases tree node inside browse tree")));
    await session.step(46, "Then the \"demo-loaded\" custom event should have fired with path \"Data Access | Databases\"", () => customFiredWith(page, "demo-loaded", "path", "Data Access | Databases"));
    await session.step(47, "And the \"Databases\" view should be current", () => viewIsCurrent(page, "Databases"));
    await session.step(48, "And the current view should be a databases view", () => currentViewType(page, "databases"));
    await session.step(49, "And no errors should have been logged", () => noErrors(page));
    await session.step(50, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Table Linking filters the demographics by the category made current", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(58, "Given Apps---Demo---Data-Access tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Data-Access tree node inside browse tree")));
    await session.step(59, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(60, "When user clicks on Apps---Demo---Data-Access---Table-Linking tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Data-Access---Table-Linking tree node inside browse tree")));
    await session.step(61, "Then the \"demo-loaded\" custom event should have fired with path \"Data Access | Table Linking\"", () => customFiredWith(page, "demo-loaded", "path", "Data Access | Table Linking"));
    await session.step(62, "And table \"Categories\" should have 8 rows", () => tableRows(page, "Categories", 8));
    await session.step(63, "And table \"Demographics\" should have 5850 rows", () => tableRows(page, "Demographics", 5850));
    await session.step(64, "When user makes row 1 of table \"Categories\" current", () => makeRowCurrentInTable(page, 1, "Categories"));
    await session.step(65, "Then 37 rows of table \"Demographics\" should pass the filter", () => tableFilterCount(page, 37, "Demographics"));
    await session.step(66, "When user makes row 7 of table \"Categories\" current", () => makeRowCurrentInTable(page, 7, "Categories"));
    await session.step(67, "Then 2444 rows of table \"Demographics\" should pass the filter", () => tableFilterCount(page, 2444, "Demographics"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
    await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Domain Databases demo turns the beta setting on, and the feature puts it back", {tag: ["@apps", "@demos", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the \"Tutorials\" package is installed", () => packageInstalled(page, "Tutorials"));
    await session.step(18, "And Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(19, "And Apps---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo tree node inside browse tree")));
    await session.step(72, "Given the \"enableDomainDatabases\" shell setting is put back at feature end", () => shellSettingPutBack(page, "enableDomainDatabases"));
    await session.step(73, "And Apps---Demo---Data-Access tree node inside browse tree is expanded", () => isExpanded(page, el("Apps---Demo---Data-Access tree node inside browse tree")));
    await session.step(74, "And user listens for \"demo-loaded\" custom event", () => listenCustom(page, "demo-loaded"));
    await session.step(75, "When user clicks on Apps---Demo---Data-Access---Domain-Databases tree node inside browse tree", () => clickOn(page, el("Apps---Demo---Data-Access---Domain-Databases tree node inside browse tree")));
    await session.step(76, "Then the \"demo-loaded\" custom event should have fired with path \"Data Access | Domain Databases\"", () => customFiredWith(page, "demo-loaded", "path", "Data Access | Domain Databases"));
    await session.step(77, "And the \"Domain Databases\" view should be current", () => viewIsCurrent(page, "Domain Databases"));
    await session.step(78, "When user clicks on \"START\" button", () => clickOn(page, el("\"START\" button")));
    await session.step(79, "Then the \"enableDomainDatabases\" shell setting should be true", () => shellSettingIs(page, "enableDomainDatabases", "true"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
});
