--name: PackageTest
--connection: System:Datagrok
select
e.friendly_name,
e.name as entity_name,
e.namespace,
e.is_deleted,
t.name as type_name,
pp.name as version_name,
p.name as package_name,
pp.version,
pp.debug
 from entities e
inner join entities_types t on t.id = e.entity_type_id
left join published_packages pp inner join packages p on p.id = pp.package_id on e.package_id = pp.id
where namespace = 'Test:' or (e.name ='Test' and t.name in ('GrokPublishedPackage', 'GrokPackage', 'Project'))
order by t.name, pp.version, e.name
