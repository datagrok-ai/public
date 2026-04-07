from typing import List, Optional
from urllib.parse import quote
from datagrok_api.http_client import HttpClient
from datagrok_api.models.package import Package, PublishedPackage


class PackagesClient:
    """
    Client for listing and retrieving Datagrok packages.

    Examples
    --------
    >>> packages = grok.packages.list()
    >>> pkg = grok.packages.get("Chem")
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def get(self, package: str) -> Package:
        """Get a package by ID or name."""
        endpoint = f"/public/v1/packages/{quote(package)}"
        response = self.client.get(endpoint)
        return Package.from_dict(response.json())

    def list(self, smart_filter: Optional[str] = None) -> List[Package]:
        """List packages with optional text filter."""
        params = {}
        if smart_filter:
            params["text"] = smart_filter
        endpoint = "/public/v1/packages"
        response = self.client.get(endpoint, params=params)
        return [Package.from_dict(data) for data in response.json()]


class PublishedPackagesClient:
    """
    Client for listing and retrieving published package versions.

    Examples
    --------
    >>> published = grok.published_packages.list()
    >>> pkg = grok.published_packages.get(published[0].id)
    """
    def __init__(self, client: HttpClient):
        self.client = client

    def get(self, id: str) -> PublishedPackage:
        """Get a published package version by ID."""
        endpoint = f"/public/v1/packages/published/{quote(id)}"
        response = self.client.get(endpoint)
        return PublishedPackage.from_dict(response.json())

    def list(self, smart_filter: Optional[str] = None,
             package_id: Optional[str] = None) -> List[PublishedPackage]:
        """List published package versions.

        Parameters
        ----------
        smart_filter : str, optional
            Text filter to search published packages.
        package_id : str, optional
            Filter by parent package ID.
        """
        params = {}
        if smart_filter:
            params["text"] = smart_filter
        if package_id:
            params["packageId"] = package_id
        endpoint = "/public/v1/packages/published"
        response = self.client.get(endpoint, params=params)
        return [PublishedPackage.from_dict(data) for data in response.json()]
