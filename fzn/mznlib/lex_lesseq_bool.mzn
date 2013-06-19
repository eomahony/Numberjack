predicate lex_lesseq_bool(array[int] of var bool: x, array[int] of var bool: y);
predicate lex_lesseq(array[int] of var bool: x, array[int] of var bool: y) = lex_lesseq_bool(x, y);
