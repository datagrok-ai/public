CREATE TABLE sal_emp (
                         name            text,
                         pay_by_quarter  integer[],
                         schedule        text[][]
);

INSERT INTO sal_emp
VALUES ('Bill',
        '{10000, 10000, 10000, 10000}',
        '{{"meeting", "lunch"}, {"training", "presentation"}}');

INSERT INTO sal_emp
VALUES ('Carol',
        '{20000, 25000, 25000, 25000}',
        '{{"breakfast", "consulting"}, {"meeting", "lunch"}}');