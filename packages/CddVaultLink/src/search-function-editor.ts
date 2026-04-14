import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ANY_READOUT_CHOICE, ANY_RUN_CHOICE, CDD_SEARCH_TYPES, CDDVaultSearchType, IN_NOT_IN_COND_CHOICES, PROTOCOL_COND_CHOICES, PROTOCOL_RUN_COND_CHOICES, ProtocolCond, ProtocolRunCond } from './constants';

const LAST_SEARCH_KEY = 'CDDVaultLink.lastSearch';

type StoredSearch = {
    protocolInNotIn?: string;
    protocolName?: string;
    runChoice?: string;
    runName?: string;
    structure?: string;
    searchType?: string;
    similarityThreshold?: number;
};
import { MoleculeFieldSearch, Protocol, queryProtocols } from './cdd-vault-api';

type CDDSerachParams = {
    structure?: string;
    structure_search_type?: CDDVaultSearchType;
    structure_similarity_threshold?: number;
    protocol?: number;
    run?: number;
}

export class SeachEditor {

    vaultId;
    protocols: { [key: string]: Protocol } = {};
    currentProtocolRuns: { [key: string]: number } = {};
    protocolNames: string[] = [];

    protocolInNotIn: DG.InputBase = ui.input.choice('Search', { value: IN_NOT_IN_COND_CHOICES[0], items: IN_NOT_IN_COND_CHOICES });
    // protocolCond: DG.InputBase = ui.input.choice('', { value: PROTOCOL_COND_CHOICES[0], items: PROTOCOL_COND_CHOICES, onValueChanged: () => {
    //     this.specificProtocolDiv.style.display = this.protocolCond.value === ProtocolCond.PROTOCOL ? 'flex' : 'none';
    // } });

    specificProtocolDiv = ui.divV([]);
    protocolListChoice: DG.InputBase = ui.input.choice('', { value: '', items: [''] });
    protocolListChoiceDiv = ui.div('', 'cdd-protocol-list-input');

    runAndReadoutDiv = ui.divV([], { style: { display: 'none' } });

    runChoice: DG.InputBase = ui.input.choice('', {
        value: PROTOCOL_RUN_COND_CHOICES[0], items: PROTOCOL_RUN_COND_CHOICES, onValueChanged: () => {
            this.runListChoiceDiv.style.display = this.runChoice.value === ProtocolRunCond.SPECIFIC ? 'flex' : 'none';
        }
    });
    runListChoiceDiv = ui.div();
    runsListChoice: DG.InputBase = ui.input.choice('', { value: '', items: [''] });

   // readoutChoice: DG.InputBase = ui.input.choice('', { value: ANY_READOUT_CHOICE, items: [ANY_READOUT_CHOICE] });
   // readoutChoiceDiv = ui.div();
   // readoutCondition: DG.InputBase = ui.input.string('');

    structureInput = ui.input.molecule('Structure');
    structureSearchTypeChoice = ui.input.choice('Search', {
        value: CDD_SEARCH_TYPES[0], items: CDD_SEARCH_TYPES, onValueChanged: () => {
            this.similarityThreshold.root.style.display = this.structureSearchTypeChoice.value === CDDVaultSearchType.SIMILARITY ? 'flex' : 'none';
        }
    });
    similarityThreshold = ui.input.float('', { value: 0.7, min: 0, max: 1, step: 0.05, showSlider: true });
    moleculeFieldsDiv = ui.divH([]);
    moleculesFields: MoleculeFieldSearch[] = [];
    plusIcon = ui.icons.add(() => { }, 'Add');

    accordion: DG.Accordion;

    initComplete: Promise<void>;
    hasRestoredSearch = false;

    constructor(vaultId: number) {
        this.vaultId = vaultId;
        this.accordion = ui.accordion(`cdd_vault_search_${this.vaultId}`);
        this.accordion.root.classList.add('cdd-search-vault-tab');
        this.protocolInNotIn.classList.add('cdd-protocol-in-not-in-input');
        //this.protocolCond.classList.add('cdd-protocol-cond-input');
        this.structureInput.classList.add('cdd-molecule-input');
        this.runListChoiceDiv.append(this.runsListChoice.root);
        this.runAndReadoutDiv.append(ui.divH([
            this.runChoice.root,
            this.runListChoiceDiv
        ], {style: {gap: '20px'}}));

        // this.readoutChoiceDiv.append(this.readoutChoice.root);
        // this.runAndReadoutDiv.append(ui.divH([
        //     this.readoutChoiceDiv,
        //     this.readoutCondition.root
        // ], {style: {gap: '20px'}}));

        this.protocolListChoiceDiv.append(this.protocolListChoice.root);
        this.specificProtocolDiv.append(this.protocolListChoiceDiv);
        this.specificProtocolDiv.append(this.runAndReadoutDiv);

        this.accordion.addPane('Protocol', () => ui.divH([
            this.protocolInNotIn.root,
            //this.protocolCond.root,
            this.specificProtocolDiv,
        ]));

        this.similarityThreshold.root.style.display = 'none';
        this.accordion.addPane('Structure', () => ui.divV([
            ui.divH([this.structureSearchTypeChoice.root, this.similarityThreshold.root], {style: {gap: '20px'}}),
            this.structureInput,
        ]));

        this.initComplete = this.init();
    }

    async init() {
        try {
            const protocolsStr = await grok.functions.call('CDDVaultLink:getProtocolsAsync', {vaultId: this.vaultId, timeoutMinutes: 1});
            const protocols = protocolsStr !== '' ? JSON.parse(protocolsStr) as Protocol[] : [];
            if (protocols && protocols.length) {
                protocols.forEach((it) => this.protocols[it.name] = it);
                this.protocolNames = Object.keys(this.protocols);
                //update the list of availabale protocols in choice input
                ui.empty(this.protocolListChoiceDiv);
                this.protocolListChoice = ui.input.choice('', {
                    value: '', items: [''].concat(this.protocolNames), onValueChanged: () => {
                        this.runAndReadoutDiv.style.display = !this.protocolListChoice!.value ? 'none' : 'flex';
                        if (this.protocolListChoice!.value) {
                           // this.updateReadoutDefinition(this.protocolListChoice!.value);
                            this.runChoice.value = ANY_RUN_CHOICE;
                            this.updateRuns(this.protocolListChoice!.value);
                        }
                    }
                });
                this.protocolListChoiceDiv.append(this.protocolListChoice.root);
            }
            this.restoreLastSearch();
        } catch(e: any) {
            grok.shell.error(e?.message ?? e);
        }
    }

    private restoreLastSearch(): void {
        const raw = grok.userSettings.getValue(LAST_SEARCH_KEY, String(this.vaultId));
        if (!raw) return;
        let s: StoredSearch;
        try {
            s = JSON.parse(raw);
        } catch { return; }
        this.hasRestoredSearch = !!(s.protocolName || s.structure);

        if (s.protocolInNotIn && IN_NOT_IN_COND_CHOICES.includes(s.protocolInNotIn as any))
            this.protocolInNotIn.value = s.protocolInNotIn;

        if (s.protocolName && this.protocolNames.includes(s.protocolName)) {
            this.protocolListChoice.value = s.protocolName;
            this.runAndReadoutDiv.style.display = 'flex';
            this.updateRuns(s.protocolName);
            if (s.runChoice && PROTOCOL_RUN_COND_CHOICES.includes(s.runChoice as any)) {
                this.runChoice.value = s.runChoice;
                this.runListChoiceDiv.style.display = s.runChoice === ProtocolRunCond.SPECIFIC ? 'flex' : 'none';
                if (s.runName && s.runName in this.currentProtocolRuns)
                    this.runsListChoice.value = s.runName;
            }
        }

        if (s.searchType && CDD_SEARCH_TYPES.includes(s.searchType as CDDVaultSearchType)) {
            this.structureSearchTypeChoice.value = s.searchType as CDDVaultSearchType;
            this.similarityThreshold.root.style.display = s.searchType === CDDVaultSearchType.SIMILARITY ? 'flex' : 'none';
        }
        if (typeof s.similarityThreshold === 'number')
            this.similarityThreshold.value = s.similarityThreshold;
        if (s.structure)
            this.structureInput.value = s.structure;
    }

    public saveLastSearch(): void {
        const payload: StoredSearch = {
            protocolInNotIn: this.protocolInNotIn.value ?? undefined,
            protocolName: this.protocolListChoice.value || undefined,
            runChoice: this.runChoice.value ?? undefined,
            runName: this.runChoice.value === ProtocolRunCond.SPECIFIC ? this.runsListChoice.value : undefined,
            structure: this.structureInput.value || undefined,
            searchType: this.structureSearchTypeChoice.value ?? undefined,
            similarityThreshold: this.similarityThreshold.value ?? undefined,
        };
        grok.userSettings.add(LAST_SEARCH_KEY, String(this.vaultId), JSON.stringify(payload));
    }

    // updateReadoutDefinition(protocolName: string) {
    //     const readouts = this.protocols[protocolName].readout_definitions.map((it) => it.name);
    //     ui.empty(this.readoutChoiceDiv);
    //     this.readoutChoice = ui.input.choice('', {
    //         value: ANY_READOUT_CHOICE, items: [ANY_READOUT_CHOICE].concat(readouts), onValueChanged: () => {
    //             this.readoutCondition.root.style.display = this.readoutChoice.value === ANY_READOUT_CHOICE ? 'none' : 'flex';
    //         }
    //     });
    //     this.readoutCondition.root.style.display = 'none';
    //     this.readoutChoiceDiv.append(this.readoutChoice.root);
    // }

    updateRuns(protocolName: string) {
        this.currentProtocolRuns = {};
        this.protocols[protocolName].runs.forEach((it) => this.currentProtocolRuns[`${it.run_date} (${it.person})`] = it.id );
        const runsNames = Object.keys(this.currentProtocolRuns);
        ui.empty(this.runListChoiceDiv);
        this.runsListChoice = ui.input.choice('', { value: runsNames.length ? runsNames[0] : '', items: runsNames});
        this.runListChoiceDiv.append(this.runsListChoice.root);
        this.runListChoiceDiv.style.display = 'none';
    }


    public getEditor(): HTMLElement {
        return this.accordion.root;
    }

    public getParams(): CDDSerachParams {
        return {
            structure: this.structureInput.value,
            structure_search_type: this.structureSearchTypeChoice.value ?? undefined,
            structure_similarity_threshold: this.similarityThreshold.value ?? undefined,
            protocol: this.protocolListChoice.value ? this.protocols[this.protocolListChoice.value].id : undefined,
            run: this.runChoice.value === ProtocolRunCond.SPECIFIC ? this.currentProtocolRuns[this.runsListChoice.value] : undefined,
        };
    }

}