select
c.name as ClientName,
ii.name as IndutryName,
at.name as AromaType,
af.name as AromaFamily,
cs.oilNettoResidualMg as oilResidual

from client c
inner join industry ii on ii.id = c.IndustryId
inner join installment i on i.ClientId = c.id
inner join service s on s.InstallmentId = i.id
inner join cartridgetoservice cs on cs.ServiceId = s.id
inner join aromatype at on at.id = cs.AromaTypeId
inner join aromafamily af on af.id = at.AromaFamilyId