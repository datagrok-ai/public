import type * as _grok from 'datagrok-api/grok';
import type * as _DG from 'datagrok-api/dg';
declare let grok: typeof _grok, DG: typeof _DG;

import {category, expect, test} from '@datagrok-libraries/test/src/test';

const layout = `{
    "#type": "ViewLayout",
    "viewStateMap": {
        "containerType": "horizontal",
        "state": {
            "width": 962,
            "height": 1185
        },
        "children": [
            {
                "containerType": "fill",
                "state": {
                    "width": 480,
                    "height": 1185,
                    "documentManager": true
                },
                "children": [
                    {
                        "containerType": "panel",
                        "state": {
                            "width": 480,
                            "height": 1185,
                            "element": {
                                "id": "31b98410-ac70-11f0-e6f2-f35d9cd08723",
                                "type": "Grid",
                                "look": {
                                    "#type": "GridLook",
                                    "showAddNewRowIcon": true,
                                    "addNewRowOnLastRowEdit": true,
                                    "topLevelDefaultMenu": true,
                                    "allowColumnMenu": true,
                                    "columns": [
                                        {
                                            "width": 28,
                                            "backgroundColor": 4294440951,
                                            "customName": "",
                                            "cellType": "row header",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 50,
                                            "columnName": "subj",
                                            "cellType": "number",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 69,
                                            "columnName": "study",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 91,
                                            "columnName": "site",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 53,
                                            "columnName": "age",
                                            "cellType": "number",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 45,
                                            "columnName": "sex",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 81,
                                            "columnName": "race",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 84,
                                            "columnName": "disease",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 89,
                                            "columnName": "started",
                                            "cellType": "number",
                                            "format": "M/d/yyyy",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 62,
                                            "columnName": "height",
                                            "cellType": "number",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 64,
                                            "columnName": "weight",
                                            "cellType": "number",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        },
                                        {
                                            "width": 66,
                                            "columnName": "control",
                                            "cellType": "bool",
                                            "headerCellStyle": {
                                                "textWrap": "words"
                                            }
                                        }
                                    ]
                                },
                                "title": "demog 10",
                                "default": true
                            }
                        },
                        "children": []
                    }
                ]
            },
            {
                "containerType": "panel",
                "state": {
                    "width": 481,
                    "height": 1185,
                    "element": {
                        "id": "31b9ab20-ac70-11f0-9495-a73b30f66fef",
                        "type": "Pie chart",
                        "look": {
                            "#type": "PieChartLook",
                            "categoryColumnName": "site",
                            "segmentAngleColumnName": "subj",
                            "legendPosition": "Auto"
                        },
                        "title": "Pie chart"
                    }
                },
                "children": []
            }
        ],
        "floating": [],
        "tableId": "2dff3b30-ac70-11f0-cec3-21abc9024bbb",
        "tableName": "demog 10"
    },
    "type": "TableView",
    "columns": [
        {
            "name": "site",
            "type": "string",
            "id": "310c5420-ac70-11f0-e8da-21d153c202d4",
            "friendlyName": "Site",
            "entityTags": [
                {
                    "tag": ".semantic-detection-duration",
                    "id": "31b9d230-ac70-11f0-de39-dbd73ee3c6c9"
                },
                {
                    "tag": ".id",
                    "id": "31b9d230-ac70-11f0-b909-ad4bd978bd72"
                }
            ],
            "metaParams": {
                ".semantic-detection-duration": "3",
                ".id": "310c5420-ac70-11f0-e8da-21d153c202d4"
            }
        },
        {
            "name": "subj",
            "type": "int",
            "id": "310c0600-ac70-11f0-df8c-c5095c75b0dd",
            "friendlyName": "Subj",
            "entityTags": [
                {
                    "tag": ".semantic-detection-duration",
                    "id": "31b9d230-ac70-11f0-ca95-b3b7bd54b277"
                },
                {
                    "tag": ".id",
                    "id": "31b9d230-ac70-11f0-9319-d7a5dd5d7c65"
                }
            ],
            "metaParams": {
                ".semantic-detection-duration": "3",
                ".id": "310c0600-ac70-11f0-df8c-c5095c75b0dd"
            }
        }
    ],
    "id": "31b9ab20-ac70-11f0-98f5-cf777ca2e738",
    "name": "Demog10",
    "friendlyName": "demog 10",
    "author": {
        "id": "878c42b0-9a50-11e6-c537-6bf8e9ab02ee"
    },
    "createdOn": "2025-11-11T11:33:07.411Z",
    "pictureId": "31b9d230-ac70-11f0-f6de-dd2e3777e1c5"
}`;

category('Dapi: layouts', () => {
  test('get applicable', async () => {
    let layout: _DG.ViewLayout| undefined;
    try {
      layout = await createTestLayout();
      const layouts = await grok.dapi.layouts.getApplicable(grok.data.demo.demog(10));
      expect(layouts.length >= 0, true, 'error in Dapi: layouts - get applicable');
    } finally {
      await safeDeleteLayout(layout);
    }

  }, { stressTest: true, owner: 'aparamonov@datagrok.ai' });

  test('filter', async () => {
    let l: _DG.ViewLayout| undefined;
    try {
      l = await createTestLayout();
      const layout = (await grok.dapi.layouts.getApplicable(grok.data.demo.demog(10)))[0];
      const layouts = (await grok.dapi.layouts.filter(`friendlyName = "${layout.friendlyName}"`).list());
      expect(layouts.length >= 0, true);
    } finally {
      await safeDeleteLayout(l);
    }

  }, { owner: 'aparamonov@datagrok.ai' });
});

async function createTestLayout(): Promise<_DG.ViewLayout> {
  const l = DG.ViewLayout.fromJson(layout);
  l.newId();
  await grok.dapi.layouts.save(l);
  return l;
}

async function safeDeleteLayout(layout?: _DG.ViewLayout): Promise<void> {
  try {
    if (layout)
      await grok.dapi.layouts.delete(layout);
  } catch (_) {}
}
