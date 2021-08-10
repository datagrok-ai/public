--name: NewUsersSummaryByDate
--input: string days
--connection: System:Datagrok
SELECT
    DATE(joined) as day, 
	COUNT(*)
FROM users
WHERE DATE(joined) BETWEEN (current_timestamp -  (CONCAT(@days, ' day'))::interval)
AND (current_timestamp - interval '0 day')
GROUP BY DATE(joined)
ORDER BY day;