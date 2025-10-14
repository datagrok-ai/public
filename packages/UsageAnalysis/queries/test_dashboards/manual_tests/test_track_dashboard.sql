--name: Manual Tests
--friendlyName: UA | Tests | Manual Tests
--connection: System:Datagrok
--input: int lastBatchesNum = 5

WITH recent_batches AS (
  SELECT DISTINCT batch_name, MAX(date_time) as latest_batch_run
  FROM test_runs r
  JOIN tests t ON t.name = r.test_name
  WHERE t.type = 'manual'
  GROUP BY batch_name
  ORDER BY latest_batch_run DESC
  LIMIT @lastBatchesNum
), recent_batches_indexed AS (
  SELECT 
  batch_name,
  ROW_NUMBER() OVER (ORDER BY latest_batch_run DESC) AS batch_index,
  latest_batch_run
  FROM recent_batches
), batch_names as (
  select batch_name from (Select distinct on ((r.params::json->>'batchName'))   
     (r.params::json->>'batchName') as batch_name,
     (r.params::json->>'start')as start,
     r.date_time as date
   from tests t full join builds b on 1 = 1
   left join test_runs r on r.test_name = t.name and r.build_name = b.name   
   where t.type = 'manual' and not (r.params::json->>'batchName') = '' 
   order by (r.params::json->>'batchName'), date) as a
   order by a.date desc
)
SELECT 
  rb.batch_name,
  rb.batch_index,
  t.name as test,
  t.type,
  r.build_name as build,
  r.date_time,
  CASE 
    WHEN r.passed IS NULL THEN 'did not run'
    WHEN r.skipped THEN 'skipped'
    WHEN r.passed THEN 'passed'
    WHEN NOT r.passed THEN 'failed'
    ELSE 'unknown'
  END as status,
  r.result
FROM tests t 
FULL JOIN recent_batches_indexed rb ON 1=1
LEFT JOIN test_runs r ON r.test_name = t.name 
  AND r.batch_name = rb.batch_name
  AND r.date_time = (
    SELECT MAX(_r.date_time) 
    FROM test_runs _r 
    WHERE _r.test_name = r.test_name 
    AND _r.batch_name = r.batch_name
  )
  LEFT JOIN (
    SELECT batch_name,
           ROW_NUMBER() OVER (ORDER BY (SELECT NULL)) AS sort_order
    FROM batch_names
) b ON r.batch_name = b.batch_name
WHERE 1=1
  AND t.type = 'manual'
  AND r.passed IS NOT NULL
  AND t.name NOT LIKE '%Unhandled exceptions: Exception'
  AND (NOT t.name = ': : ')
ORDER  BY b.sort_order, r.batch_name, t.name
--end