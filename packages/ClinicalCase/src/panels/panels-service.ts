import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import { CLIN_TRIAL_GOV_SEARCH, HttpService } from '../services/http.service';

export async function createPropertyPanel(viewClass: any) {
    switch (viewClass.name) {
        case 'Summary': {
            grok.shell.o = await summaryPanel(viewClass.studyId);
            break;
        }
        default: {
            break;
        }
    }
}

export async function summaryPanel(studyId: string) {
    const httpService = new HttpService();
    let clinTrialsGovInfo = await httpService.getStudyData('R01NS050536');
    //let clinTrialsGovInfo = await httpService.getStudyData(studyId);
    let studyLink = `${CLIN_TRIAL_GOV_SEARCH}${clinTrialsGovInfo['NCTId']}`
    clinTrialsGovInfo[`Study link`] = ui.link('Go to study page', ()=>{window.open(studyLink, '_blank').focus();})
    let root = ui.div();
    root.appendChild(ui.h2(`Study ${studyId} summary`));
    root.appendChild(ui.tableFromMap(clinTrialsGovInfo));
    return root;
}