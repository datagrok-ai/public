import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { studyConfigJsonFileName } from '../constants/constants';
import { SUMMARY_VIEW_NAME } from '../constants/view-names-constants';
import { StudyConfig } from '../types/types';
import { studyLoadedSubject, studies, openStudy, createStudyWithConfig, addStudyToBrowseTree } from '../utils/app-utils';
import { studyConfigToMap } from '../utils/utils';

export function createImportStudyView(): DG.ViewBase {
  const importView = DG.View.create();
  importView.root.classList.add('preclinical-case-study-import-view');
  importView.name = 'Import Study - Preclinical Case';

  let studyConfig: StudyConfig | null = null;
  const clinicalCaseNode = grok.shell.browsePanel.mainTree.getOrCreateGroup('Apps').
    getOrCreateGroup('Preclinical Case');

  const fileNamesDiv = ui.div();
  const errorDiv = ui.div();
  const filesSection = ui.divV([
    fileNamesDiv,
    errorDiv,
  ], {style: {marginTop: '10px'}});
  const statisticsDiv = ui.div('', {style: {marginLeft: '10px', marginTop: '10px', position: 'relative'}});
  const importStudyStatusDiv = ui.div('', {style: {position: 'relative', marginLeft: '20px'}});

  const filesInput = ui.input.files('', {
    onValueChanged: async () => {
      ui.empty(fileNamesDiv);
      ui.empty(errorDiv);
      ui.empty(statisticsDiv);
      ui.empty(importStudyStatusDiv);

      if (filesInput.value && filesInput.value.length > 0) {
        const fileTagsContainer = ui.div(
          filesInput.value.map((file) => {
            const tag = ui.div(file.name, {
              style: {
                display: 'inline-block',
                padding: '4px 8px',
                marginRight: '8px',
                marginBottom: '4px',
                backgroundColor: '#E7F0F3',
                color: 'var(--blue-6)',
                borderRadius: '4px',
                fontSize: '12px',
              },
            });
            tag.classList.add('preclinical-case-import-study-file-tag');
            return tag;
          }),
          {style: {display: 'flex', flexWrap: 'wrap', marginTop: '10px', marginLeft: '20px'}},
        );
        fileNamesDiv.append(fileTagsContainer);

        try {
          ui.setUpdateIndicator(statisticsDiv, true, 'Collecting study summary...');
          importStudyButton.disabled = true;

          studyConfig = await createStudyWithConfig(filesInput.value, clinicalCaseNode, true);

          if (studyConfig) {
            const configMap = studyConfigToMap(studyConfig);
            const statsTable = ui.tableFromMap(configMap);
            statisticsDiv.append(ui.divText('Summary', {style: {marginLeft: '10px'}}));
            statisticsDiv.append(statsTable);
            ui.setUpdateIndicator(statisticsDiv, false);
            importStudyButton.disabled = false;
          }
        } catch (e: any) {
          const errorMessage = e?.message ?? String(e);
          errorDiv.append(ui.divText(`Error: ${errorMessage}`,
            {style: {color: 'red', marginLeft: '20px', marginTop: '10px'}}));
          ui.empty(statisticsDiv);
          importStudyButton.disabled = true;
          grok.shell.error(e);
        }
      }
    },
    acceptExtensions: ['.xpt', '.csv', '.xml', '.json'],
  });
  filesInput.root.style.paddingTop = '6px';

  const importStudyButton = ui.button('Import', async () => {
    importStudyStatusDiv.classList.add('preclinical-case-study-import-in-progress-div');
    ui.setUpdateIndicator(importStudyStatusDiv, true, `Loading data for study ${studyConfig!.name}`);
    const sub = studyLoadedSubject.subscribe((data) => {
      if (data.name === studyConfig!.name) {
        sub.unsubscribe();
        if (data.loaded) {
          importStudyButton.disabled = true;
          importStudyStatusDiv.append(ui.divText(`Study ${studyConfig!.name} loaded successfully`,
            {style: {color: 'green'}}));
          const tags = Array.from(fileNamesDiv.querySelectorAll('.preclinical-case-import-study-file-tag'));
          if (data.errorDomains) {
            for (const errorDomain of data.errorDomains) {
              const tag = tags.filter((it) => it.textContent.toLocaleLowerCase() === errorDomain.toLowerCase());
              if (tag.length)
                (tag[0] as HTMLElement).style.color = 'red';
            }
          }
        } else {
          importStudyStatusDiv.append(ui.divText(`Error loading study ${studyConfig!.name}`,
            {style: {color: 'red'}}));
        }
        importStudyStatusDiv.classList.remove('preclinical-case-study-import-in-progress-div');
        ui.setUpdateIndicator(importStudyStatusDiv, false);
      }
    });
    for (const file of filesInput.value!) {
      await grok.dapi.files.write(
        `System:AppData/Preclinicalcase/SEND/${studyConfig!.name}/${file.name}`,
        Array.from(file.data));
    }
    grok.dapi.files.writeAsText(
      `System:AppData/Preclinicalcase/SEND/${studyConfig!.name}/${studyConfigJsonFileName}`,
      JSON.stringify(studyConfig));
    addStudyToBrowseTree(studies[studyConfig!.name!], clinicalCaseNode, filesInput.value!);
    openStudy(clinicalCaseNode, studyConfig!.name!, SUMMARY_VIEW_NAME);
  });

  importStudyButton.disabled = true;

  importView.root.append(ui.divV([
    ui.divH([ui.h3('Import study files'),
      filesInput.root, importStudyButton], {style: {gap: '10px', marginLeft: '20px'}}),
    importStudyStatusDiv,
    filesSection,
    statisticsDiv,
  ]));

  return importView;
}