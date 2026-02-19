//name: TestFunctionScriptJoinWithoutCols
//language: javascript
//input: dataframe data [demog.csv]
//output: dataframe result {action:join(data); meta: {"StudyTag": {"description": "Arbitrary study-level tag/label for demo purposes", "Tag.Category": "Demo", "Tag.Source": "Script"}, "SiteCode": {"description": "Arbitrary site code for demo purposes", "Tag.Category": "Demo", "Tag.Source": "Script"}, "QCFlag": {"description": "Arbitrary QC flag for demo purposes", "Tag.Category": "Demo", "Tag.Source": "Script"}, "Notes": {"description": "Free-text notes for demo purposes", "Tag.Category": "Demo", "Tag.Source": "Script"}}} [Columns to join back to the original table]

result = DG.DataFrame.create(data.rowCount);

const studyTag = DG.Column.string('StudyTag', data.rowCount).init((i) => `STUDY-${(i % 3) + 1}`);
const siteCode = DG.Column.int('SiteCode', data.rowCount).init((i) => 100 + (i % 20));
const qcFlag = DG.Column.bool('QCFlag', data.rowCount).init((i) => (i % 7) !== 0);
const notes = DG.Column.string('Notes', data.rowCount).init((i) => (i % 10 === 0) ? 'Check source record' : '');

studyTag.setTag('description', 'Arbitrary study-level tag/label for demo purposes');
studyTag.setTag('Tag.Category', 'Demo');
studyTag.setTag('Tag.Source', 'Script');

siteCode.setTag('description', 'Arbitrary site code for demo purposes');
siteCode.setTag('Tag.Category', 'Demo');
siteCode.setTag('Tag.Source', 'Script');

qcFlag.setTag('description', 'Arbitrary QC flag for demo purposes');
qcFlag.setTag('Tag.Category', 'Demo');
qcFlag.setTag('Tag.Source', 'Script');

notes.setTag('description', 'Free-text notes for demo purposes');
notes.setTag('Tag.Category', 'Demo');
notes.setTag('Tag.Source', 'Script');

result.columns.add(studyTag);
result.columns.add(siteCode);
result.columns.add(qcFlag);
result.columns.add(notes);
