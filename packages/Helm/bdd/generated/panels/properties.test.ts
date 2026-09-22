/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/panels/properties.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [helm.panel.properties]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {helmInitialized, setLongPeptide} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {isExpanded, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnUnits, currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Properties pane of a HELM cell", () => {
  const session = feature(test, "features/panels/properties.feature", import.meta.url);
  test("The Properties pane of a HELM cell", {tag: ["@journey", "@realizes:helm.panel.properties"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the Helm package is initialized", () => helmInitialized(page));
    await session.step(13, "And user opens helm-showcase dataset", () => openDataset(page, ds("helm-showcase")));
    await session.step(14, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(15, "Then \"HELM\" column should have units \"helm\"", () => columnUnits(page, "HELM", "helm"));
    await run.scenario("The current cell's formula, weight and extinction coefficient", async () => {
      await session.step(18, "When user clicks on the \"cell 1 of HELM\" area of grid", () => clickArea(page, "cell 1 of HELM", el("grid")));
      await session.step(19, "Then row 1 should be current", () => currentRowIs(page, 1));
      await session.step(20, "And \"Properties\" section in context panel should be visible", () => shouldBe(page, el("\"Properties\" section in context panel"), "visible"));
      await session.step(21, "Given \"Properties\" section in context panel is expanded", () => isExpanded(page, el("\"Properties\" section in context panel")));
      await session.step(22, "Then \"Properties\" section in context panel should contain text \"C6H12N2O3S\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "C6H12N2O3S"));
      await session.step(23, "And \"Properties\" section in context panel should contain text \"192.23\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "192.23"));
      await session.step(24, "And \"Properties\" section in context panel should contain text \"0.06\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "0.06"));
      await session.step(25, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Another cell shows its own values", async () => {
      await session.step(29, "When user clicks on the \"cell 2 of HELM\" area of grid", () => clickArea(page, "cell 2 of HELM", el("grid")));
      await session.step(30, "Then row 2 should be current", () => currentRowIs(page, 2));
      await session.step(31, "And \"Properties\" section in context panel should contain text \"C50H77N13O15S\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "C50H77N13O15S"));
      await session.step(32, "And \"Properties\" section in context panel should contain text \"1132.30\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "1132.30"));
      await session.step(33, "And \"Properties\" section in context panel should not contain text \"C6H12N2O3S\"", () => shouldNotContainText(page, el("\"Properties\" section in context panel"), "C6H12N2O3S"));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A sequence over 1000 characters gets a warning, not a calculation", async () => {
      await session.step(37, "When user sets \"HELM\" column in row 5 to a peptide of 600 alanines", () => setLongPeptide(page, "HELM", 5, 600));
      await session.step(38, "And user clicks on the \"cell 5 of HELM\" area of grid", () => clickArea(page, "cell 5 of HELM", el("grid")));
      await session.step(39, "Then row 5 should be current", () => currentRowIs(page, 5));
      await session.step(40, "And \"Properties\" section in context panel should contain text \"Too long sequence\"", () => shouldContainText(page, el("\"Properties\" section in context panel"), "Too long sequence"));
      await session.step(41, "And \"Properties\" section in context panel should not contain text \"formula\"", () => shouldNotContainText(page, el("\"Properties\" section in context panel"), "formula"));
      await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
