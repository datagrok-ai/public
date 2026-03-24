from datagrok_api.models import Package, PublishedPackage


def test_list_packages(grok):
    packages = grok.packages.list()
    assert isinstance(packages, list)
    assert all(isinstance(p, Package) for p in packages)


def test_get_package_by_id(grok):
    packages = grok.packages.list()
    if packages:
        pkg = grok.packages.get(packages[0].id)
        assert pkg.id == packages[0].id
        assert pkg.name is not None


def test_get_package_by_name(grok):
    packages = grok.packages.list()
    if packages:
        pkg = grok.packages.get(packages[0].name)
        assert pkg.name == packages[0].name


def test_list_packages_with_filter(grok):
    packages = grok.packages.list()
    if packages:
        name = packages[0].name
        filtered = grok.packages.list(smart_filter=f'shortName = "{name}"')
        assert len(filtered) >= 1
        assert any(p.name == name for p in filtered)


def test_list_published_packages(grok):
    published = grok.published_packages.list()
    assert isinstance(published, list)
    assert all(isinstance(p, PublishedPackage) for p in published)


def test_get_published_package(grok):
    published = grok.published_packages.list()
    if published:
        pkg = grok.published_packages.get(published[0].id)
        assert pkg.id == published[0].id
        assert pkg.version is not None


def test_list_published_by_package_id(grok):
    packages = grok.packages.list()
    if packages:
        published = grok.published_packages.list(package_id=packages[0].id)
        assert isinstance(published, list)
        assert all(isinstance(p, PublishedPackage) for p in published)
