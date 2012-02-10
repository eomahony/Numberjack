var x0 integer, >=0, <= 9;
var x1, >=0, <= 6;

param V > 0;

maximize objective:
    V * x0 + 15000 * x1;

subject to st0:
    0.3 * x0 + 0.4 * x1 >= 2.0;

subject to st1:
    0.4 * x0 + 0.2 * x1 >= 1.5;

subject to st2:
    0.2 * x0 + 0.3 * x1 >= 0.5;
