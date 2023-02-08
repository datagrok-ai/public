--name: NewEventsSummaryByDate
--meta.cache: true
--input: string days
--connection: System:Datagrok
SELECT
    DATE(event_time) as day, 
	COUNT(*)
FROM events
WHERE DATE(event_time) BETWEEN (current_timestamp -  (CONCAT(@days, ' day'))::interval)
AND (current_timestamp - interval '0 day')
GROUP BY DATE(event_time)
ORDER BY day;