// @ts-check
/** @type {import('@docusaurus/plugin-content-docs').SidebarsConfig} */
const typedocSidebar = {
  "items": [
    {
      "type": "category",
      "label": "datagrok_client",
      "link": {
        "type": "doc",
        "id": "py/datagrok_client/index"
      },
      "items": [
        {
          "type": "category",
          "label": "Classes",
          "items": [
            {
              "type": "doc",
              "id": "py/datagrok_client/classes/DatagrokClient",
              "label": "DatagrokClient"
            }
          ]
        }
      ]
    },
    {
      "type": "category",
      "label": "group",
      "link": {
        "type": "doc",
        "id": "py/group/index"
      },
      "items": [
        {
          "type": "category",
          "label": "Classes",
          "items": [
            {
              "type": "doc",
              "id": "py/group/classes/Group",
              "label": "Group"
            },
            {
              "type": "doc",
              "id": "py/group/classes/GroupRelation",
              "label": "GroupRelation"
            },
            {
              "type": "doc",
              "id": "py/group/classes/GroupMembershipRequest",
              "label": "GroupMembershipRequest"
            }
          ]
        }
      ]
    },
    {
      "type": "category",
      "label": "model",
      "link": {
        "type": "doc",
        "id": "py/model/index"
      },
      "items": [
        {
          "type": "category",
          "label": "Classes",
          "items": [
            {
              "type": "doc",
              "id": "py/model/classes/Model",
              "label": "Model"
            }
          ]
        }
      ]
    }
  ]
};
module.exports = typedocSidebar.items;
