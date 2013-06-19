predicate lex_less_bool(array[int] of var bool: x, array[int] of var bool: y);
predicate lex_less(array[int] of var bool: x, array[int] of var bool: y) = lex_less_bool(x, y);
