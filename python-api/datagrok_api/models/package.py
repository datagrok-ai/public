from typing import Optional
from datetime import datetime
from datagrok_api.models.model import NamedModel


class Package(NamedModel):
    """
    Represents a Datagrok package (plugin).

    Attributes
    ----------
    full_name : Optional[str]
        Full package name (e.g., npm scope name).
    description : Optional[str]
        Package description.
    category : Optional[str]
        Package category.
    image_url : Optional[str]
        URL of the package icon/image.
    url : Optional[str]
        Package repository or homepage URL.
    """
    def __init__(self,
                 name: str,
                 id: Optional[str] = None,
                 friendly_name: Optional[str] = None,
                 full_name: Optional[str] = None,
                 description: Optional[str] = None,
                 category: Optional[str] = None,
                 image_url: Optional[str] = None,
                 url: Optional[str] = None,
                 created_on: Optional[datetime] = None,
                 updated_on: Optional[datetime] = None):
        super().__init__(id=id, name=name, friendly_name=friendly_name,
                         created_on=created_on, updated_on=updated_on)
        self.full_name = full_name
        self.description = description
        self.category = category
        self.image_url = image_url
        self.url = url

    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "friendlyName": self.friendly_name,
            "fullName": self.full_name,
            "description": self.description,
            "category": self.category,
            "imageUrl": self.image_url,
            "url": self.url,
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'Package':
        pkg = cls(
            name=data.get("name", ""),
            id=data.get("id"),
            friendly_name=data.get("friendlyName"),
            full_name=data.get("fullName"),
            description=data.get("description"),
            category=data.get("category"),
            image_url=data.get("imageUrl"),
            url=data.get("url"),
        )
        created_str = data.get("createdOn")
        pkg.created_on = datetime.fromisoformat(created_str) if created_str else None
        updated_str = data.get("updatedOn")
        pkg.updated_on = datetime.fromisoformat(updated_str) if updated_str else None
        return pkg


class PublishedPackage(NamedModel):
    """
    Represents a published version of a Datagrok package.

    Attributes
    ----------
    version : Optional[str]
        Version string (e.g., "1.2.3").
    is_current : bool
        Whether this is the currently active version.
    debug : bool
        Whether this is a debug (dev) build.
    package_author : Optional[str]
        Author of the package.
    published_on : Optional[datetime]
        When this version was published.
    commit : Optional[str]
        Git commit hash.
    """
    def __init__(self,
                 name: str,
                 id: Optional[str] = None,
                 friendly_name: Optional[str] = None,
                 version: Optional[str] = None,
                 is_current: bool = False,
                 debug: bool = False,
                 package_author: Optional[str] = None,
                 published_on: Optional[datetime] = None,
                 commit: Optional[str] = None,
                 created_on: Optional[datetime] = None,
                 updated_on: Optional[datetime] = None):
        super().__init__(id=id, name=name, friendly_name=friendly_name,
                         created_on=created_on, updated_on=updated_on)
        self.version = version
        self.is_current = is_current
        self.debug = debug
        self.package_author = package_author
        self.published_on = published_on
        self.commit = commit

    def to_dict(self):
        return {
            "id": self.id,
            "name": self.name,
            "friendlyName": self.friendly_name,
            "version": self.version,
            "isCurrent": self.is_current,
            "debug": self.debug,
            "packageAuthor": self.package_author,
            "publishedOn": self.published_on.isoformat() if self.published_on else None,
            "commit": self.commit,
            "createdOn": self.created_on.isoformat() if self.created_on else None,
            "updatedOn": self.updated_on.isoformat() if self.updated_on else None,
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'PublishedPackage':
        pkg = cls(
            name=data.get("name", ""),
            id=data.get("id"),
            friendly_name=data.get("friendlyName"),
            version=data.get("version"),
            is_current=data.get("isCurrent", False),
            debug=data.get("debug", False),
            package_author=data.get("packageAuthor"),
            commit=data.get("commit"),
        )
        for field, key in [("created_on", "createdOn"), ("updated_on", "updatedOn"), ("published_on", "publishedOn")]:
            val = data.get(key)
            if val:
                setattr(pkg, field, datetime.fromisoformat(val))
        return pkg
