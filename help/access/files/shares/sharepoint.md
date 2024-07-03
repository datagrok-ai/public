---
title: "SharePoint"
---

Provides access to the
[SharePoint](https://www.microsoft.com/en-us/microsoft-365/sharepoint) or [OneDrive](https://www.microsoft.com/en-us/microsoft-365/onedrive) as
[file shares](../files.md).

## Application settings

Datagrok connects with SharePoint through Microsoft Entra ID, that needs a registered application in Entra ID. Application client ID and application secret should be copied from Entra ID to **Settings/Server/Connectors**.

## Configuration

Connector needs to specify **Domain**, which is an URL of a sharepoint instance. **Site** specifies which site under sharepoint domain is used. Usually it is a **root site**, in every other case site url's explicitly specifies name in the format of `https://test.sharepoint.com/sites/test` instead of a domain `https://iekonstantinamelichev.sharepoint.com` as-is. Drive is a name of a document library that will be used as a file share.


Sharepoint connector uses OAuth for authorization. When configuring connection, click on **"connect..."** link near **Access token** field.

![Sharepoint configuration](../img/sharepoint-parameters.png)

## OneDrive

SharePoint connector supports OneDrive as well, when **Drive** is set to "OneDrive".


## Connection parameters

````json
{
  "parameters": {
    "domain": "",
    "root site": true,
    "site": "", // specify if not a root site
    "drive": "",
    "redirect url": "" // is set automatically
  },
  "credentials": {
    "parameters": {
      "access token": "", // is set automatically
      "refresh token": "", // is set automatically
    }
  }
}
````

See also:

* [File shares](../files.md)
* [Microsoft Entra ID](https://learn.microsoft.com/en-us/azure/databricks/dev-tools/app-aad-token)
* [SharePoint](https://www.microsoft.com/en-us/microsoft-365/sharepoint)
* [OneDrive](https://www.microsoft.com/en-us/microsoft-365/onedrive)
