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

    get steps() { return 18; }

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
        'Analyze Peptides', this.getMenuItem('Bio', true));
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

        await this.action('Click the Gear Icon (⚙)', gearIconClickEvent, gearIcon as HTMLElement);
        
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
        await this.action('Set Similarity Threshold to 90', similarityChangedPromise, similarityThresholdInput as HTMLElement);

        const webGpuCheckbox: HTMLInputElement | null = d.root.querySelector('input[name="input-Use-WebGPU"]');

        if (webGpuCheckbox && !webGpuCheckbox.disabled) {
            const webGpuChangedPromise = new Promise<void>((resolve) => {
                webGpuCheckbox?.addEventListener('click', () => {
                    if (webGpuCheckbox.checked)
                        resolve();
                })
            })
            await this.action('Check Use WebGPU', webGpuChangedPromise, webGpuCheckbox as HTMLElement);
        }

        okBtn.disabled = false;
        // click ok button on the dialog
        await this.action('Click OK to start analysis', d.onClose, okBtn);
        
        // wait for SAR view to be created
        await this.action('Wait for analysis to complete',
              grok.events.onViewerAdded.pipe(operators.filter((data: DG.EventData) => {
                const found = data.args.viewer.type === 'Logo Summary Table';
                return found;
              })));
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
                this._removeHints(step3Hint);
                step3Hint.remove();
                resolve();
            });
            step3NextButton.addEventListener('click', () => {
                this._removeHints(step3Hint);
                step3Hint.remove();
                resolve();
            });
        });
        await this.action('Click NEXT to proceed', step3Prom);

        // ########## Step 4: Monomer grid ##########
        this.title('Monomer Grid', true);
        this.describe(`The peptide column is split into positions; each cell is a monomer in that position.<br>
             Hover a monomer to preview its chemical structure.`);
        const step4Hint = greenHint(gridRoot, paragraphs([`Hover over a <i>monomer cell</i> in the <i>main grid</i> to see its structure`]),'right');
        await this.action('Hover over a monomer cell in the main table grid',
            grok.events.onTooltipShown.pipe(operators.filter(() => {
                return ui.tooltip?.root?.textContent?.toLowerCase()?.includes('helmcore') ?? false;
            })));

        this._removeHints(step4Hint);
        step4Hint.remove();



        // ########## Step 5: WebLogo (select a monomer@position) ##########
        this.title('WebLogo Header', true);
        this.describe('Above each position, the WebLogo shows which monomers occur and how often. Click a letter to select that monomer@position across all viewers.');
    
        const step5Hint = greenHint(gridRoot, paragraphs(['Hover over WebLogo header to preview statistics.','Click a letter to select that <b>monomer@position</b>.']), 'top');
        const webLogoClick =
            new Promise<void>((res) => {
                const timer = setInterval(() => {
                    if (tv.dataFrame.selection.trueCount > 0) {
                        clearInterval(timer);
                        res();
                    }
                }, 200);
                setTimeout(() => {
                    clearInterval(timer);
                    res();
                }, 1800000); // 30 minutes timeout
            })

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
        const step6B1Prom = new Promise<void>((resolve) => {
          okBtn1.addEventListener('click', () => resolve());
            attachRemovingHintListener(step6B1, () => resolve());
        });
        await this.action('Click NEXT to proceed', step6B1Prom);
        this._removeHints(step6B1);
        step6B1.remove();

        // Balloon 2 — SAR view root, press Esc to clear selection
        const step6B2 = greenHint(tv.root, paragraphs(['All viewers are synchronized and show your <i>WebLogo</i> selection.',' Press <b>Esc</b> to clear the selection.']), 'top');
        const escPromise =  new Promise<void>((resolve) => {
          const onKey = (e: KeyboardEvent) => {
            if (e.key === 'Escape') {
              const sub = tv.dataFrame.selection.onChanged.subscribe(() => {
                if (tv.dataFrame.selection.trueCount === 0) {
                  sub.unsubscribe();
                  document.removeEventListener('keydown', onKey);
                  resolve();
                }
              });
              if (tv.dataFrame.selection.trueCount === 0) {
                sub.unsubscribe();
                document.removeEventListener('keydown', onKey);
                resolve();
              }
            }
          };
          document.addEventListener('keydown', onKey);
        });
        await this.action('Press Esc to clear selection', escPromise);
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
        const svmGearClick: Promise<void> = (svmGear ? rxjs.fromEvent(svmGear, 'click') : rxjs.fromEvent(svmRoot, 'click')).pipe(operators.take(1), operators.map(() => void 0)).toPromise();
        await this.action('Open SVM settings (gear)', svmGearClick, svmGear ?? svmRoot);
        this._removeHints(step7B1);
        step7B1.remove();

        // Balloon 2 — context panel shows SVM settings
        const nextBtn = ui.button('NEXT', () => {});
        const step7B2Content = ui.divV([
          paragraphs(['The <i>Context Panel</i> contains settings for the <i>Sequence Variability Map</i>.']),
          ui.divH([nextBtn])
        ]);
        const step7B2 = greenHint(contextPanelRoot, step7B2Content, 'left');
        const step7b2Prom = new Promise<void>((resolve) => {
            nextBtn.addEventListener('click', () => resolve());
            attachRemovingHintListener(step7B2, () => resolve());
        });
        await this.action('Click NEXT to proceed', step7b2Prom);
        this._removeHints(step7B2);
        step7B2.remove();

        // ########## Step 8: Sequence Variability Map · Mutation Cliffs ##########
        const s81b = paragraphs(['Mutation Cliffs', 'In this mode, each <b>cell</b> shows <b>counts</b> of sequence pairs that differ only at that <br>monomer(row)-position(column) (<b>Size</b>) and <br>their mean activity difference (<b>Color</b>).','<b>Hover/Click</b> any non-empty cell.']);
        const step8B1 = greenHint(svmRoot, s81b, 'top');
        const mutViewerClickPromise = new Promise<void>((res) => {
            const timer = setInterval(() => {
                const o = grok.shell.o;
                if (o instanceof HTMLElement && o.getElementsByClassName('d4-accordion-title').length > 0 &&
                Array.from(o.getElementsByClassName('d4-accordion-title')).some((el) => el.textContent?.toLowerCase()?.includes('selection sources') &&
                el.textContent?.toLowerCase()?.includes('mutation cliffs')) && o.getElementsByClassName('d4-pane-mutation_cliffs_pairs').length > 0)
                {
                    clearInterval(timer);
                    res();
                }
            }, 200);
            setTimeout(() => {
                clearInterval(timer);
                res();
            }, 1800000); // 30 minutes timeout
        });
        await this.action('Click a Mutation Cliffs cell', mutViewerClickPromise
        );
        this._removeHints(step8B1);
        step8B1.remove();

        
        const mutationCliffsPanel = (grok.shell.o as HTMLElement).getElementsByClassName('d4-pane-mutation_cliffs_pairs')[0] as HTMLElement;
        if (mutationCliffsPanel) {
            const step8B2 = greenHint(mutationCliffsPanel, paragraphs(['Mutation Cliffs context panel shows sequence pairs only differring at selected position, their activity distributions, and more']), 'left');
            const step8B2Ok = ui.button('NEXT', () => {});
            step8B2.appendChild(step8B2Ok);
            const mutPanelNextPromise = new Promise<void>((resolve) => {step8B2Ok.addEventListener('click', () => resolve())
                attachRemovingHintListener(step8B2, () => resolve());
            });
            await this.action('Click NEXT to proceed', mutPanelNextPromise);

            this._removeHints(step8B2);
            step8B2.remove();
        }
        
        // ########## Step 9: Sequence Variability Map · Invariant map ##########
        const invariantMapInputRoot = Array.from(svmRoot.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelectorAll('.ui-input-bool.ui-input-root') ?? [])
            .find((r) => (r.textContent ?? '').toLowerCase().includes('invariant map'))!;
        const invariantRadioButton = invariantMapInputRoot.querySelector('input[type="radio"]') as HTMLInputElement;
        const invariantPromise = new Promise<void>((res) => {
            const int = setInterval(() => {
                if ((invariantRadioButton as HTMLInputElement).checked == true || invariantRadioButton.value == 'true') {
                    clearInterval(int);
                    res();
                }
            }, 200);
            setTimeout(() => {
                clearInterval(int);
                res();
            }, 1800000); // 30 minutes timeout
        });
        const invariantHint = greenHint(invariantRadioButton, paragraphs(['Switch the <i>SVM viewer</i> to <i>Invariant Map</i> mode using the radio button.']), 'top');
        await this.action('Switch SVM mode to Invariant Map', invariantPromise, invariantRadioButton);
        this._removeHints(invariantHint);
        invariantHint.remove();

        // Step 9: Invariant Map
        const step91Hint = greenHint(svmRoot, paragraphs(['Invariant Map', 'In this mode, each <b>cell</b> shows how many sequences contain <br>that monomer at that position (<b>number</b>) and <br>the mean activity of those sequences (<b>color</b>).','<b>Hover</b> or <b>click</b> any cell to see details']), 'right');
        const step91Next = ui.button('NEXT', () => {});
        step91Hint.appendChild(step91Next);
        const step91Prom = new Promise<void>((resolve) => {
            step91Next.addEventListener('click', () => resolve());
            attachRemovingHintListener(step91Hint, () => resolve());
        });
        await this.action('Click NEXT to proceed', step91Prom);
        this._removeHints(step91Hint);
        step91Hint.remove();

        // ########## Step 10: Most Potent Residues (orientation) ##########
        const mpr = Array.from(tv.viewers).find((v: DG.Viewer) => v.type === 'Most Potent Residues')!;

        const mprRoot = (mpr ? mpr.root : tv.root);
        this.title('Most Potent Residues', true);
        this.describe('The <b>Most Potent Residues</b> viewer highlights the most potent monomers at each position, based on statistical distributions. It helps identify key monomers contributing to high/low activity. Use the Gear icon to adjust settings.');

        const step10Hint = greenHint(mprRoot, paragraphs(['The <i>Most Potent Residues</i> viewer highlights the most potent monomers at each position.','Use the <b>gear</b> icon to adjust its settings.']), 'left');
        const mprGear = mprRoot.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;
        const mprGearClick: Promise<void> = (mprGear ? rxjs.fromEvent(mprGear, 'click') : rxjs.fromEvent(mprRoot, 'click')).pipe(operators.take(1), operators.map(() => void 0)).toPromise();
        await this.action('Open Most Potent Residues settings (gear)', mprGearClick, mprGear ?? mprRoot);
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
        const mclNextProm = new Promise<void>((resolve) => {
            step10MCLNext.addEventListener('click', () => resolve());
            attachRemovingHintListener(step10MCLHint, () => resolve());
        });
        await this.action('Click NEXT to proceed', mclNextProm);
        this._removeHints(step10MCLHint);
        step10MCLHint.remove();

        const step10LSTHint = greenHint(LSTViewer.root, paragraphs(['The Logo Summary Table details <i>clusters</i>, generates their <i>WebLogos</i>, along with other statistics.',' Click the <b>gear</b> icon to adjust its settings.']), 'left');
        const lstGear = LSTViewer.root.parentElement!.parentElement!.parentElement!.parentElement!.parentElement!.querySelector('.grok-icon.grok-font-icon-settings') as HTMLElement;

        const lstGearClick: Promise<void> = (lstGear ? rxjs.fromEvent(lstGear, 'click') : rxjs.fromEvent(LSTViewer.root, 'click')).pipe(operators.take(1), operators.map(() => void 0)).toPromise();
        const lstContextPromise = new Promise<void>(async (resolve) => {
            await lstGearClick;
            const int = setInterval(() => {
                if (grok.shell.o === LSTViewer) {
                    clearInterval(int);
                    resolve();
                }
        }, 200);
        setTimeout(() => {
            clearInterval(int);
            resolve();
        }, 1800000); // 30 minutes timeout
    });
        await this.action('Open Logo Summary Table settings (gear)', lstContextPromise, lstGear ?? LSTViewer.root);
        this._removeHints(step10LSTHint);
        step10LSTHint.remove();

        // expand all sections in LST settings
        (Array.from(contextPanelRoot.getElementsByClassName('property-grid-category-body-hide') ?? []) as HTMLElement[]).forEach((el) => el.classList.remove('property-grid-category-body-hide'));

        const aggColumnsHost = Array.from(contextPanelRoot.getElementsByClassName('property-grid-multi-column-editor') ?? []).find((el) => el.textContent?.toLowerCase()?.includes('0 / 25'));

        const aggColumnsButton = aggColumnsHost!.querySelector('button') as HTMLButtonElement;
        const aggColumnsClick = rxjs.fromEvent(aggColumnsButton, 'click').pipe(operators.take(1)).toPromise();
        const lstHint = greenHint(contextPanelRoot, paragraphs(['Add aggregated column','Under <b>Aggregations</b>, choose position <b>14</b> (Column named "14") and hit <b>OK</b> to add a <b>pie chart</b> distribution column.']), 'left');
        const lstAggRegationPromise = new Promise<void>(async (resolve) => {
            await aggColumnsClick;
            const int = setInterval(() => {
                if (LSTViewer.props.columns.length > 0) {
                    clearInterval(int);
                    resolve();
                }
            }, 200);
            setTimeout(() => {
                clearInterval(int);
                resolve();
            }, 1800000); // 30 minutes timeout
        });

        await this.action('Add pie chart aggregation for position 14', lstAggRegationPromise, aggColumnsButton);
        this._removeHints(lstHint);
        lstHint.remove();

        // scroll through the added columns to see the pie chart
        const lstHorzScroll = LSTViewer.root.querySelector('.d4-range-selector.d4-grid-horz-scroll') as HTMLElement;
        if (lstHorzScroll) {

            const mouseRelease = rxjs.fromEvent(lstHorzScroll, 'mousedown').pipe(operators.take(1),  operators.map(() => void 0)).toPromise();
            const scrollHint = greenHint(lstHorzScroll, paragraphs(['Scroll the Logo Summary Table horizontally to see the newly added pie chart column.']), 'bottom');
            const nextBtn = ui.button('NEXT', () => {});
            scrollHint.appendChild(nextBtn);
            const scrollPromise = new Promise<void>((resolve) => {
                mouseRelease.then(() => resolve());
                nextBtn.addEventListener('click', () => resolve());
                attachRemovingHintListener(scrollHint, () => resolve());
            });
            await this.action('Scroll horizontally in Logo Summary Table. Click NEXT to proceed to next step', scrollPromise, lstHorzScroll);
            this._removeHints(scrollHint);
            scrollHint.remove();
        }

        this.title('Update Analysis Configuration', true);
        this.describe(`You can update the analysis configuration by clicking the Wrench icon in the top-right corner of the SAR view.<br>
            This allows you to adjust parameters and re-run the analysis to see how changes impact the results.`);
        
        const pepAnalWrench = document.querySelector('.fal.fa-wrench') as HTMLElement;
        // const pepAnalWrenchClick = rxjs.fromEvent(pepAnalWrench, 'click').pipe(operators.take(1), operators.map(() => void 0)).toPromise();
        const wrenchHint = greenHint(pepAnalWrench, paragraphs(['Lastly, you can update your SAR settings and analysis anytime.','Click the <b>wrench</b> icon, enable <b>Dendrogram</b>, then <b>OK</b> to add it to the analysis.']), 'bottom');
        const analDialog = await this.openDialog('Click the Wrench to update analysis configuration', 'Peptides settings', pepAnalWrench);
        this._removeHints(wrenchHint);
        wrenchHint.remove();
        // expand all sections in settings
        analDialog.root.querySelectorAll('.d4-accordion-pane-content').forEach((el) => (el as HTMLElement).classList.add('expanded'));
        const dendrogramCheckBox = analDialog.root.querySelector('input[name="input-Dendrogram"]') as HTMLInputElement;
        const dendrogramPromise = new Promise<void>((res) => {
            const int = setInterval(() => {
                if (dendrogramCheckBox.checked) {
                    clearInterval(int);
                    res();
                }
            }, 200);
            setTimeout(() => {
                clearInterval(int);
                res();
            }, 1800000); // 30 minutes timeout
        });
        await this.action('Check Dendrogram', dendrogramPromise, dendrogramCheckBox);
        const analDialogOk = analDialog.root.querySelector('button.ui-btn.ui-btn-ok') as HTMLButtonElement;
        await this.action('Click OK to re-run analysis', analDialog.onClose, analDialogOk);

        // wait for dendrogram to appear for 2 seconds
        await new Promise<void>((resolve) => {
            setTimeout(() => resolve(), 2000);
        });

        const finalHintDiv = ui.divV([paragraphs([`All set!`, `You launched <b>SAR</b>, explored <b>WebLogos</b> and <b>monomers</b>,<br>configured <b>Sequence Variability Map</b> modes,<br>reviewed <b>Most Potent Residues</b>, examined <b>clusters</b>,<br>enriched the <b>Logo Summary Table</b>, and updated <b>global settings</b>.`]),
            ui.link('Read more', 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar')
        ])
        const finalHint = greenHint(svmRoot, finalHintDiv, 'top');
        const finalOk = ui.button('OK', () => {});
        finalHint.appendChild(finalOk);
        await new Promise<void>((resolve) => {
            finalOk.addEventListener('click', () => resolve());
            attachRemovingHintListener(finalHint, () => resolve());
        });
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

function paragraphs(texts: string[]) {
    return ui.divV(texts.map((t, i) => {
        const thing = i == 0 ? ui.h2(t, {style: {marginBottom: '0.8em'}}) : ui.divText(t);
        thing.innerHTML = t;
        return thing;
}));
}

function attachRemovingHintListener(hintDiv: HTMLElement, onRemoved: () => void) {
    const timer = setInterval(() => {
        if (!document.body.contains(hintDiv)) {
            clearInterval(timer);
            onRemoved();
        }
    }, 200);
    setTimeout(() => {clearInterval(timer); onRemoved();}, 1800000); // 30 miunutes timeout
}