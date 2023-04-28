CREATE TABLE operators (id bigint, json_data jsonb, path_data path, circle_data circle);

INSERT INTO operators (id, json_data, path_data, circle_data) VALUES (1, '{"reading": 1.230e-5}', path('((0,0),(1,0))'), circle('((0,0),10)'));

INSERT INTO operators (id, json_data, path_data, circle_data) VALUES (2, '{"bar": "baz", "balance": 7.77, "active":false}', path(
        '(15.878137629895164,47.08306448089695),
         (15.56169808311181,47.219041634920686),
         (15.267442604782124,47.4201665137259),
         (15.092631384557304,47.71366328136526),
         (15.234428926980286,47.95865145177352)'
    ), circle('((0,0),15)'));
