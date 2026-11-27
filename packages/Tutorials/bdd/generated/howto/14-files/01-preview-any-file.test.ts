/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/14-files/01-preview-any-file.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelClosed, contextPanelClosed, packageInstalled, simpleModeOff, toolboxPaneHidden, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Files", () => {
  const session = feature(test, "features/howto/14-files/01-preview-any-file.feature", import.meta.url);
  test("Preview any file from Browse", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And the \"FileEditors\" package is installed", () => packageInstalled(page, "FileEditors"));
    await session.step(18, "And the \"Chem\" package is installed", () => packageInstalled(page, "Chem"));
    await session.step(19, "And the \"BiostructureViewer\" package is installed", () => packageInstalled(page, "BiostructureViewer"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(22, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(23, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "When the video says \"Browse > Files > Demo: sample files of many kinds\"", () => videoSays(page, "Browse > Files > Demo: sample files of many kinds"));
    await session.step(25, "And user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(26, "And user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(27, "And user expands Files---Demo tree node inside browse tree", () => expand(page, el("Files---Demo tree node inside browse tree")));
    await session.step(28, "And user clicks on Files---Demo---cars.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(29, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(30, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(31, "When the video says \"One click shows a file here, without opening it\"", () => videoSays(page, "One click shows a file here, without opening it"));
    await session.step(32, "And user expands Files---Demo---images tree node inside browse tree", () => expand(page, el("Files---Demo---images tree node inside browse tree")));
    await session.step(33, "And user expands Files---Demo---images---cats tree node inside browse tree", () => expand(page, el("Files---Demo---images---cats tree node inside browse tree")));
    await session.step(34, "And user clicks on Files---Demo---images---cats---cat2.jpg tree node inside browse tree", () => clickOn(page, el("Files---Demo---images---cats---cat2.jpg tree node inside browse tree")));
    await session.step(35, "Then the \"cat2.jpg\" view should be current", () => viewIsCurrent(page, "cat2.jpg"));
    await session.step(36, "When the video says \"Images show as they are\"", () => videoSays(page, "Images show as they are"));
    await session.step(37, "And user collapses Files---Demo---images tree node inside browse tree", () => collapse(page, el("Files---Demo---images tree node inside browse tree")));
    await session.step(38, "And user expands Files---Demo---texts tree node inside browse tree", () => expand(page, el("Files---Demo---texts tree node inside browse tree")));
    await session.step(39, "And user expands Files---Demo---texts---pdf tree node inside browse tree", () => expand(page, el("Files---Demo---texts---pdf tree node inside browse tree")));
    await session.step(40, "And user clicks on Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree", () => clickOn(page, el("Files---Demo---texts---pdf---research-paper.pdf tree node inside browse tree")));
    await session.step(41, "Then the \"research-paper.pdf\" view should be current", () => viewIsCurrent(page, "research-paper.pdf"));
    await session.step(42, "When the video says \"A PDF shows its pages\"", () => videoSays(page, "A PDF shows its pages"));
    await session.step(43, "And user collapses Files---Demo---texts tree node inside browse tree", () => collapse(page, el("Files---Demo---texts tree node inside browse tree")));
    await session.step(44, "And user expands Files---Demo---chem tree node inside browse tree", () => expand(page, el("Files---Demo---chem tree node inside browse tree")));
    await session.step(45, "And user expands Files---Demo---chem---sdf tree node inside browse tree", () => expand(page, el("Files---Demo---chem---sdf tree node inside browse tree")));
    await session.step(46, "And user clicks on Files---Demo---chem---sdf---ApprovedDrugs2015.sdf tree node inside browse tree", () => clickOn(page, el("Files---Demo---chem---sdf---ApprovedDrugs2015.sdf tree node inside browse tree")));
    await session.step(47, "Then the \"ApprovedDrugs2015\" view should be current", () => viewIsCurrent(page, "ApprovedDrugs2015"));
    await session.step(48, "When the video says \"A molecule file (SDF) shows its structures\"", () => videoSays(page, "A molecule file (SDF) shows its structures"));
    await session.step(49, "And user collapses Files---Demo---chem tree node inside browse tree", () => collapse(page, el("Files---Demo---chem tree node inside browse tree")));
    await session.step(50, "And user expands Files---Demo---bio tree node inside browse tree", () => expand(page, el("Files---Demo---bio tree node inside browse tree")));
    await session.step(51, "And user expands Files---Demo---bio---ngl-formats tree node inside browse tree", () => expand(page, el("Files---Demo---bio---ngl-formats tree node inside browse tree")));
    await session.step(52, "And user clicks on Files---Demo---bio---ngl-formats---1crn.pdb tree node inside browse tree", () => clickOn(page, el("Files---Demo---bio---ngl-formats---1crn.pdb tree node inside browse tree")));
    await session.step(53, "Then the \"1crn.pdb\" view should be current", () => viewIsCurrent(page, "1crn.pdb"));
    await session.step(54, "And \"Reset Camera\" button should be visible", () => shouldBe(page, el("\"Reset Camera\" button"), "visible"));
    await session.step(55, "And no table should be left in the workspace", () => noTableLeft(page));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
    await session.step(57, "When the video says \"A protein structure (PDB) shows in 3D\"", () => videoSays(page, "A protein structure (PDB) shows in 3D"));
  });
});
