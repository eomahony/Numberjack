predicate table_bool(array[int] of var bool: x, array[int, int] of bool: t);
predicate table(array[int] of var bool: x, array[int, int] of bool: t) = table_bool(x, t);

predicate table_bool_reif(array[int] of var bool: x, array[int, int] of bool: t, var bool: b) = abort("Reified table/2 for Booleans is not supported.");
