predicate lex_less_int(array[int] of var int: x, array[int] of var int: y);
predicate lex_less(array[int] of var int: x, array[int] of var int: y) = lex_less_int(x,y);
