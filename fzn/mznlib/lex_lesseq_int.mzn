predicate lex_lesseq_int(array[int] of var int: x, array[int] of var int: y);
predicate lex_lesseq(array[int] of var int: x, array[int] of var int: y) = lex_lesseq_int(x, y);
