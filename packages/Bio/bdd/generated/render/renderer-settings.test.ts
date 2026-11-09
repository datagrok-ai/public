/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/render/renderer-settings.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [bio.panel.sequence-renderer, bio.rendering.biln]
--- */
import {test} from '@playwright/test';
import '../../bindings/elements.js';
import '../../bindings/monomer-form.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {bioInitialized} from '../../bindings/steps.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, expand, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, columnUnits} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {contextPanelOpen, contextPanelShows, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sequence renderer settings on the context panel", () => {
  const session = feature(test, "features/render/renderer-settings.feature", import.meta.url);
  test("Show Connections starts checked on a BILN column and toggling it writes the column's tag", {tag: ["@realizes:bio.panel.sequence-renderer", "@realizes:bio.rendering.biln"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the Bio package is initialized", () => bioInitialized(page));
    await session.step(19, "And user opens BILN dataset", () => openDataset(page, ds("BILN")));
    await session.step(20, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(21, "Then \"biln\" column should have units \"biln\"", () => columnUnits(page, "biln", "biln"));
    await session.step(24, "When user clicks on the \"header biln\" area of grid", () => clickArea(page, "header biln", el("grid")));
    await session.step(25, "Then the context panel should show \"biln\"", () => contextPanelShows(page, "biln"));
    await session.step(26, "When user expands Bioinformatics accordion header in context panel", () => expand(page, el("Bioinformatics accordion header in context panel")));
    await session.step(27, "And user expands \"Sequence Renderer\" accordion header in context panel", () => expand(page, el("\"Sequence Renderer\" accordion header in context panel")));
    await session.step(28, "And user expands \"Renderer Settings\" accordion header in context panel", () => expand(page, el("\"Renderer Settings\" accordion header in context panel")));
    await session.step(29, "Then \"Show Connections\" checkbox in context panel should be checked", () => shouldBe(page, el("\"Show Connections\" checkbox in context panel"), "checked"));
    await session.step(30, "When user unchecks \"Show Connections\" checkbox in context panel", () => uncheck(page, el("\"Show Connections\" checkbox in context panel")));
    await session.step(31, "Then \"biln\" column should have tag \".mm.cellRenderer.showConnections\" equal to \"false\"", () => columnTag(page, "biln", ".mm.cellRenderer.showConnections", "false"));
    await session.step(32, "When user checks \"Show Connections\" checkbox in context panel", () => check(page, el("\"Show Connections\" checkbox in context panel")));
    await session.step(33, "Then \"biln\" column should have tag \".mm.cellRenderer.showConnections\" equal to \"true\"", () => columnTag(page, "biln", ".mm.cellRenderer.showConnections", "true"));
    await session.step(34, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(35, "And no errors should have been logged", () => noErrors(page));
  });
});
