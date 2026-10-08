import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Tutorial, TutorialPrerequisites } from "@datagrok-libraries/tutorials/src/tutorial";
import { getPlatform, Platform } from "../../shortcuts";
import { _package } from '../../../package';
import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import * as operators from 'rxjs/operators';

export class PeptidesSarTutorial extends Tutorial {

    get name() {
        return 'Peptides SAR';
    }

    get description() {
        return `This tutorial demonstrates how to analyze structure-activity relationships (SAR) in peptide datasets.
         We will explore peptide data, visualize SAR patterns based on monomer-positions, cluster data based on biological distances,
         identify mutation cliffs and compute statistical distributions of activity values.`;
    }

    // 24 as launched here; a machine with WebGPU lists one more (Check Use WebGPU), which the counter clamps
    get steps() { return 24; }

    get icon() {
        return '🧬📈';
    }

    helpUrl: string = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    prerequisites: TutorialPrerequisites = {packages: ['Bio', 'Peptides']};
    platform: Platform = getPlatform();

    protected async _run() {
        grok.shell.windows.showBrowse = false;
        grok.shell.windows.showContextPanel = false;
        grok.shell.windows.showProjects = false;
        grok.shell.windows.showToolbox = false;
        this.header.textContent = this.name;
        //grok.shell.windows.showContextPanel = true;
        this.describe(`The <b>Peptides SAR</b> tool detects and visualizes pairs of peptides
            with highly similar sequence but significantly different activity levels.
            The resulting view consists of various viewers, which detail monomer-position based statistical activity distributions,
            clusters and mutation cliffs<hr>`);
        this.t = await grok.data.loadTable(`${_package.webRoot}files/MSA.csv`);
        const tv = grok.shell.addTableView(this.t);
        this.title('Launch SAR analysis', true);

        this.describe(`When you open a biological dataset, Datagrok automatically detects sequences, along with their notation (e.g., FASTA, Separator, HELM etc.),
            and shows sequence-specific tools, actions, and information. Access them through:<br>
            <ul>
            <li><b>Bio</b> menu (it contains all Bioinformatics tools)</li>
            </ul><br>
            Let's launch the Peptides SAR tool.`);
        
        const d = await this.openDialog('On the Top Menu, click Bio > Analyze > SAR...',
        'Analyze Peptides', this.getMenuItem('Bio', true), '',
        () => grok.shell.v.ribbonMenu.find('Bio | Analyze | SAR...').click());
        const okBtn = $(d.root).find('button.ui-btn.ui-btn-ok')[0] as HTMLButtonElement;
        okBtn.disabled = true;

        const gearIcon = d.root.querySelector('.grok-icon.grok-font-icon-settings');
        if (!gearIcon)
            throw new Error('Cannot find settings icon in SAR dialog');
        const gearIconClickEvent = (async () => {await rxjs.fromEvent(gearIcon, 'click').pipe(operators.take(1)).toPromise();})();
        
        // ########## Step 2: Adjust clustering settings ##########
        this.title('Adjust Clustering Settings', true);

        this.describe(`In the <b>Peptides SAR</b> dialog, you can specify the sequence column, activity column, scaling method and clustering parameters.<br>
            Clustering parameters include distance metric, threshold, and others. <br>
            Lets adjust the Similarity threshold to 90 and enable WebGPU acceleration for faster clustering computation.`); 

        await this.action('Click the Gear Icon (⚙)', gearIconClickEvent, gearIcon as HTMLElement, '',
            () => (gearIcon as HTMLElement).click());
        
        const similarityThresholdInput: HTMLInputElement | null = d.root.querySelector('input[name="input-Similarity-Threshold"]');
        if (!similarityThresholdInput)
            throw new Error('Cannot find Similarity Threshold input in SAR dialog');

        const similarityChangedPromise = new Promise<void>((resolve) => {
            similarityThresholdInput.addEventListener('input', () => {
                if (similarityThresholdInput.value?.toString()?.split('.')?.[0] === '90')
                    resolve();
            })
            similarityThresholdInput.addEventListener('change', () => {
                if (similarityThresholdInput.value?.toString()?.split('.')?.[0] === '90')
                    resolve();
            })
        })
        await this.action('Set Similarity Threshold to 90', similarityChangedPromise, similarityThresholdInput as HTMLElement, '',
            () => Tutorial.setInputValue(similarityThresholdInput, '90'));

        const webGpuCheckbox: HTMLInputElement | null = d.root.querySelector('input[name="input-Use-WebGPU"]');

        if (webGpuCheckbox && !webGpuCheckbox.disabled) {
            const webGpuChangedPromise = new Promise<void>((resolve) => {
                webGpuCheckbox?.addEventListener('click', () => {
                    if (webGpuCheckbox.checked)
                        resolve();
                })
            })
            await this.action('Check Use WebGPU', webGpuChangedPromise, webGpuCheckbox as HTMLElement, '',
                () => webGpuCheckbox.click());
        }

        okBtn.disabled = false;
        // click ok button on the dialog
        await this.action('Click OK to start analysis', d.onClose, okBtn, '', () => okBtn.click());
        
        // wait for SAR view to be created
        await this.action('Wait for analysis to complete',
              grok.events.onViewerAdded.pipe(operators.filter((data: DG.EventData) => {
                const found = data.args.viewer.type === 'Logo Summary Table';
                return found;
              })), null, '', null);
        grok.shell.windows.showContextPanel = true;
        

        // ########## Step 3: Explore SAR View ##########
        const grid = tv.grid;
        const gridRoot = grid.root;
        const step3Hint = greenHint(grid.root, paragraphs([`You don't see the original table, but it's still open.`,
            `All original columns remain available in <b>filters</b>/<b>selectors</b>
            alongside newly generated ones corresponding to position-wise monomers and scaled activity.`]),'right');
        const step3NextButton = ui.button('NEXT', () => {});
            step3Hint.appendChild(step3NextButton);
        const step3Hintesc = step3Hint.querySelector('.fa-times') as HTMLElement;
        this._placeHints(step3NextButton);
        const step3Prom = new Promise<void>((resolve) => {
            step3Hintesc.addEventListener('click', () => {
                resolve();
            });
            step3NextButton.addEventListener('click', () => {
                resolve();
            });
        });
        await this.action('Click NEXT to proceed', step3Prom, null, '', () => step3NextButton.click());
        this._removeHints(step3Hint);
        step3Hint.remove();

        // ########## Step 4: Monomer grid ##########
        this.title('Monomer Grid', true);
        this.describe(`The peptide column is split into positions; each cell is a monomer in that position.<br>
             Hover a monomer to preview its chemical structure.`);
        const step4Hint = greenHint(gridRoot, paragraphs([`Hover over a <i>monomer cell</i> in the <i>main grid</i> to see its structure`]),'right');
        // a monomer's tooltip names the library it comes from, which differs between stands
        await this.action('Hover over a monomer cell in the main table grid',
            grok.events.onTooltipShown.pipe(operators.filter(() =>
                gridRoot.matches(':hover') && ui.tooltip?.root?.querySelector('.ui-form.ui-tooltip') != null)));

        this._removeHints(step4Hint);
        step4Hint.remove();



        // ########## Step 5: WebLogo (select a monomer@position) ##########
        this.title('WebLogo Header', true);
        this.describe('Above each position, the WebLogo shows which monomers occur and how often. Click a letter to select that monomer@position across all viewers.');
    
        const step5Hint = greenHint(gridRoot, paragraphs(['Hover over WebLogo header to preview statistics.','Click a letter to select that <b>monomer@position</b>.']), 'top');
        const webLogoClick = poll(() => tv.dataFrame.selection.trueCount > 0);

        await this.action('Click any WebLogo letter', webLogoClick);
        this._removeHints(step5Hint);
        step5Hint.remove();

        // ########## Step 6: Explore monomers for a position (selection feedback) ##########
        // Balloon 1 — Context Panel
        grok.shell.windows.showContextPanel = true;
        const contextPanelRoot = document.querySelector('.grok-prop-panel') as HTMLElement ?? grok.shell.v.root;
        const okBtn1 = ui.button('NEXT', () => {});
        const step6B1Content = ui.divV([
          paragraphs(['Context Panel shows information and statistics for your selected <i>monomer-positions</i>.']),
          ui.divH([okBtn1])
        ]);
        const step6B1 = greenHint(contextPanelRoot, step6B1Content, 'left');
        const step6B1Prom = nextOrHintGone(okBtn1, step6B1);
        await this.action('Click NEXT to proceed', step6B1Prom, null, '', () => okBtn1.click());
        this._removeHints(step6B1);
        step6B1.remove();

        // Balloon 2 — SAR view root, press Esc to clear selection
        const step6B2 = greenHint(tv.root, paragraphs(['All viewers are synchronized and show your <i>WebLogo</i> selection.',' Press <b>Esc</b> to clear the selection.']), 'top');
        const escPromise = rxjs.fromEvent<KeyboardEvent>(document, 'keydown').pipe(
          operators.filter((e) => e.key === 'Escape'),
          operators.switchMap(() => rxjs.concat(rxjs.of(null), tv.dataFrame.selection.onChanged)),
          operators.filter(() => tv.dataFrame.selection.trueCount === 0));
        await this.action('Press Esc to clear selection', escPromise, null, '',
            Tutorial.apiSkip(() => tv.dataFrame.selection.setAll(false)));
        this._removeHints(step6B2);
        step6B2.remove();

        this.title('Sequence Variability Map', true);
        this.describe(`The <b>Sequence Variability Map</b> viewer shows mutation distributions and statistics across all positions.
            It helps identify positions where mutations significantly impact activity, guiding further analysis. <br>
            The viewer operates in two modes: Mutation Cliffs and Invariant Map. <br>`);

        // ########## Step 7: Sequence Variability Map (open settings) ##########
        const svm = Array.from(tv.viewers).find((v: DG.Viewer) => v.type === 'Sequence Variability Map');
        const svmRoot = (svm ? svm.root : tv.root);
        const step7B1 = greenHint(svmRoot, paragraphs(['The <i>Sequence Variability Map</i> shows mutation distribution.', 'Click the <b>Gear icon</b> to adjust its settings.']), 'right');
        const svmGear = svmRoot.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;
        const svmGearClick = rxjs.fromEvent(svmGear ?? svmRoot, 'click');
        await this.action('Open SVM settings (gear)', svmGearClick, svmGear ?? svmRoot, '', () => svmGear.click());
        this._removeHints(step7B1);
        step7B1.remove();

        // Balloon 2 — context panel shows SVM settings
        const nextBtn = ui.button('NEXT', () => {});
        const step7B2Content = ui.divV([
          paragraphs(['The <i>Context Panel</i> contains settings for the <i>Sequence Variability Map</i>.']),
          ui.divH([nextBtn])
        ]);
        const step7B2 = greenHint(contextPanelRoot, step7B2Content, 'left');
        const step7b2Prom = nextOrHintGone(nextBtn, step7B2);
        await this.action('Click NEXT to proceed', step7b2Prom, null, '', () => nextBtn.click());
        this._removeHints(step7B2);
        step7B2.remove();

        // ########## Step 8: Sequence Variability Map · Mutation Cliffs ##########
        const s81b = paragraphs(['In this mode, each <b>cell</b> shows <b>counts</b> of sequence pairs that differ only at that monomer(row)-position(column) (<b>Size</b>) and their mean activity difference (<b>Color</b>).','<b>Hover/Click</b> any non-empty cell.'], 'Mutation Cliffs');
        const step8B1 = greenHint(svmRoot, s81b, 'top');
        const mutViewerClickPromise = poll(() => {
            const o = grok.shell.o;
            return o instanceof HTMLElement && o.getElementsByClassName('d4-accordion-title').length > 0 &&
                Array.from(o.getElementsByClassName('d4-accordion-title')).some((el) => el.textContent?.toLowerCase()?.includes('selection sources') &&
                el.textContent?.toLowerCase()?.includes('mutation cliffs')) && o.getElementsByClassName('d4-pane-mutation_cliffs_pairs').length > 0;
        });
        await this.action('Click a Mutation Cliffs cell', mutViewerClickPromise
        );
        this._removeHints(step8B1);
        step8B1.remove();

        
        const mutationCliffsPanel = grok.shell.o instanceof HTMLElement ?
            grok.shell.o.getElementsByClassName('d4-pane-mutation_cliffs_pairs')[0] as HTMLElement : null;
        if (mutationCliffsPanel) {
            const step8B2 = greenHint(mutationCliffsPanel, paragraphs(['Mutation Cliffs context panel shows sequence pairs only differring at selected position, their activity distributions, and more']), 'left');
            const step8B2Ok = ui.button('NEXT', () => {});
            step8B2.appendChild(step8B2Ok);
            const mutPanelNextPromise = nextOrHintGone(step8B2Ok, step8B2);
            await this.action('Click NEXT to proceed', mutPanelNextPromise, null, '', () => step8B2Ok.click());

            this._removeHints(step8B2);
            step8B2.remove();
        }
        
        // ########## Step 9: Sequence Variability Map · Invariant map ##########
        const invariantMapInputRoot = Array.from(svmRoot.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelectorAll('.ui-input-bool.ui-input-root') ?? [])
            .find((r) => (r.textContent ?? '').toLowerCase().includes('invariant map'))!;
        const invariantRadioButton = invariantMapInputRoot.querySelector('input[type="radio"]') as HTMLInputElement;
        const invariantPromise = poll(() => invariantRadioButton.checked || invariantRadioButton.value == 'true');
        const invariantHint = greenHint(invariantRadioButton, paragraphs(['Switch the <i>SVM viewer</i> to <i>Invariant Map</i> mode using the radio button.']), 'top');
        await this.action('Switch SVM mode to Invariant Map', invariantPromise, invariantRadioButton, '',
            () => invariantRadioButton.click());
        this._removeHints(invariantHint);
        invariantHint.remove();

        // Step 9: Invariant Map
        const step91Hint = greenHint(svmRoot, paragraphs(['In this mode, each <b>cell</b> shows how many sequences contain that monomer at that position (<b>number</b>) and the mean activity of those sequences (<b>color</b>).','<b>Hover</b> or <b>click</b> any cell to see details'], 'Invariant Map'), 'right');
        const step91Next = ui.button('NEXT', () => {});
        step91Hint.appendChild(step91Next);
        const step91Prom = nextOrHintGone(step91Next, step91Hint);
        await this.action('Click NEXT to proceed', step91Prom, null, '', () => step91Next.click());
        this._removeHints(step91Hint);
        step91Hint.remove();

        // ########## Step 10: Most Potent Residues (orientation) ##########
        const mpr = Array.from(tv.viewers).find((v: DG.Viewer) => v.type === 'Most Potent Residues')!;

        const mprRoot = (mpr ? mpr.root : tv.root);
        this.title('Most Potent Residues', true);
        this.describe('The <b>Most Potent Residues</b> viewer highlights the most potent monomers at each position, based on statistical distributions. It helps identify key monomers contributing to high/low activity. Use the Gear icon to adjust settings.');

        const step10Hint = greenHint(mprRoot, paragraphs(['The <i>Most Potent Residues</i> viewer highlights the most potent monomers at each position.','Use the <b>gear</b> icon to adjust its settings.']), 'left');
        const mprGear = mprRoot.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;
        const mprGearClick = rxjs.fromEvent(mprGear ?? mprRoot, 'click');
        await this.action('Open Most Potent Residues settings (gear)', mprGearClick, mprGear ?? mprRoot, '',
            () => mprGear.click());
        this._removeHints(step10Hint);
        step10Hint.remove();


        // Step 10: MCL
        const mclViewer = Array.from(tv.viewers).find((v: DG.Viewer) => v.type === 'MCL')!;
        const LSTViewer = Array.from(tv.viewers).find((v: DG.Viewer) => v.type === 'Logo Summary Table')!;
        this.title('Clusters and Logo Summary Table', true);
        this.describe(`During the analysis, sequences are <b>clustered</b> based on the selected distance function and algorithm (MCL/UMAP/T-SNE).<br>
            These clusters are visualized in the <b>MCL scatterplot</b> viewer and used for generation of the <b>Logo Summary Table</b>, which details each cluster, their activity distributions along with other statistics.`);
        const step10MCLHint = greenHint(mclViewer.root, paragraphs(['This <i>scatterplot</i> shows <i>clusters</i> based on the selected <i>algorithm</i>.','To learn how it works, complete the <b>scatterplot tutorial</b>.']), 'left');
        const step10MCLNext = ui.button('NEXT', () => {});
        step10MCLHint.appendChild(step10MCLNext);
        const mclNextProm = nextOrHintGone(step10MCLNext, step10MCLHint);
        await this.action('Click NEXT to proceed', mclNextProm, null, '', () => step10MCLNext.click());
        this._removeHints(step10MCLHint);
        step10MCLHint.remove();

        const step10LSTHint = greenHint(LSTViewer.root, paragraphs(['The Logo Summary Table details <i>clusters</i>, generates their <i>WebLogos</i>, along with other statistics.',' Click the <b>gear</b> icon to adjust its settings.']), 'left');
        const lstGear = LSTViewer.root.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;

        const lstContextPromise = rxjs.fromEvent(lstGear ?? LSTViewer.root, 'click').pipe(operators.take(1),
            operators.switchMap(() => poll(() => grok.shell.o === LSTViewer)));
        await this.action('Open Logo Summary Table settings (gear)', lstContextPromise, lstGear ?? LSTViewer.root, '',
            () => (lstGear ?? LSTViewer.root).click());
        this._removeHints(step10LSTHint);
        step10LSTHint.remove();

        // expand all sections in LST settings
        (Array.from(contextPanelRoot.getElementsByClassName('property-grid-category-body-hide') ?? []) as HTMLElement[]).forEach((el) => el.classList.remove('property-grid-category-body-hide'));

        const aggColumnsHost = Array.from(contextPanelRoot.getElementsByClassName('property-grid-multi-column-editor') ?? []).find((el) => el.textContent?.toLowerCase()?.includes('0 / 25'));

        const aggColumnsButton = aggColumnsHost!.querySelector('button') as HTMLButtonElement;
        const lstHint = greenHint(contextPanelRoot, paragraphs(['Under <b>Aggregations</b>, choose position <b>14</b> (Column named "14") and hit <b>OK</b> to add a <b>pie chart</b> distribution column.'], 'Add aggregated column'), 'left');
        const lstAggRegationPromise = rxjs.fromEvent(aggColumnsButton, 'click').pipe(operators.take(1),
            operators.switchMap(() => poll(() => LSTViewer.props.columns.length > 0)));

        await this.action('Add pie chart aggregation for position 14', lstAggRegationPromise, aggColumnsButton, '',
            Tutorial.apiSkip(() => LSTViewer.setOptions({columns: ['14']})));
        this._removeHints(lstHint);
        lstHint.remove();

        // scroll through the added columns to see the pie chart
        const lstHorzScroll = LSTViewer.root.querySelector('.d4-range-selector.d4-grid-horz-scroll') as HTMLElement;
        if (lstHorzScroll) {

            const scrollHint = greenHint(lstHorzScroll, paragraphs(['Scroll the Logo Summary Table horizontally to see the newly added pie chart column.']), 'bottom');
            const nextBtn = ui.button('NEXT', () => {});
            scrollHint.appendChild(nextBtn);
            const scrollPromise = rxjs.merge(rxjs.fromEvent(lstHorzScroll, 'mousedown'), nextOrHintGone(nextBtn, scrollHint));
            await this.action('Scroll horizontally in Logo Summary Table. Click NEXT to proceed to next step',
                scrollPromise, lstHorzScroll, '', () => nextBtn.click());
            this._removeHints(scrollHint);
            scrollHint.remove();
        }

        this.title('Update Analysis Configuration', true);
        this.describe(`You can update the analysis configuration by clicking the Wrench icon in the top-right corner of the SAR view.<br>
            This allows you to adjust parameters and re-run the analysis to see how changes impact the results.`);
        
        const pepAnalWrench = document.querySelector('.fal.fa-wrench') as HTMLElement;
        // const pepAnalWrenchClick = rxjs.fromEvent(pepAnalWrench, 'click').pipe(operators.take(1), operators.map(() => void 0)).toPromise();
        const wrenchHint = greenHint(pepAnalWrench, paragraphs(['Lastly, you can update your SAR settings and analysis anytime.','Click the <b>wrench</b> icon, enable <b>Dendrogram</b>, then <b>OK</b> to add it to the analysis.']), 'bottom');
        const analDialog = await this.openDialog('Click the Wrench to update analysis configuration', 'Peptides settings', pepAnalWrench,
            '', () => pepAnalWrench.click());
        this._removeHints(wrenchHint);
        wrenchHint.remove();
        // expand all sections in settings
        analDialog.root.querySelectorAll('.d4-accordion-pane-content').forEach((el) => (el as HTMLElement).classList.add('expanded'));
        const dendrogramCheckBox = analDialog.root.querySelector('input[name="input-Dendrogram"]') as HTMLInputElement;
        const dendrogramPromise = poll(() => dendrogramCheckBox.checked);
        await this.action('Check Dendrogram', dendrogramPromise, dendrogramCheckBox, '',
            () => dendrogramCheckBox.click());
        const analDialogOk = analDialog.root.querySelector('button.ui-btn.ui-btn-ok') as HTMLButtonElement;
        await this.action('Click OK to re-run analysis', analDialog.onClose, analDialogOk, '',
            () => analDialogOk.click());

        // wait for dendrogram to appear for 2 seconds
        await new Promise<void>((resolve) => {
            setTimeout(() => resolve(), 2000);
        });

        const finalHintDiv = ui.divV([paragraphs([`You launched <b>SAR</b>, explored <b>WebLogos</b> and <b>monomers</b>, configured <b>Sequence Variability Map</b> modes, reviewed <b>Most Potent Residues</b>, examined <b>clusters</b>, enriched the <b>Logo Summary Table</b>, and updated <b>global settings</b>.`], `All set!`),
            ui.link('Read more', 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar')
        ])
        const finalHint = greenHint(svmRoot, finalHintDiv, 'top');
        const finalOk = ui.button('OK', () => {});
        finalHint.appendChild(finalOk);
        await this.firstEvent(nextOrHintGone(finalOk, finalHint));
        this._removeHints(finalHint);
        finalHint.remove();
        
        // THE END

    }
}

function greenHint(el: HTMLElement, content: HTMLElement, position?: "top" | "bottom" | "left" | "right" | undefined): HTMLElement {
    const hint = ui.hints.addHint(el, content, position);
    hint.classList.add('ui-hint-popup-green-1');
    return hint;
}

function paragraphs(texts: string[], title?: string) {
    const items: HTMLElement[] = [];
    if (title != null) {
        const heading = ui.h2(title, {style: {marginBottom: '0.8em'}});
        heading.innerHTML = title;
        items.push(heading);
    }
    for (const t of texts) {
        const p = ui.divText(t);
        p.innerHTML = t;
        items.push(p);
    }
    return ui.divV(items);
}

// observables rather than promises, since the tutorial unsubscribes a step's stream when it is closed, so no poll outlives it
function poll(condition: () => boolean): rxjs.Observable<number> {
    return rxjs.interval(200).pipe(operators.filter(() => condition()));
}

function nextOrHintGone(next: HTMLElement, hint: HTMLElement): rxjs.Observable<unknown> {
    return rxjs.merge(rxjs.fromEvent(next, 'click'), poll(() => !document.body.contains(hint)));
}