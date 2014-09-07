/*=============================================================================
 * parser for CSP instances represented in XML format
 * 
 * Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *=============================================================================
 */
#ifndef _XMLParser_constants_h_
#define _XMLParser_constants_h_

/*
 * declarations shared by C and C++
 */

/* different definitions of relations */
typedef enum {REL_CONFLICT,REL_SUPPORT,REL_SOFT} RelType;

/* different ways to represent an expression */
typedef enum {TREE,PREFIX,INFIX_C,POSTFIX,MATHML} Syntax;

/* different type of nodes in an AST (abstract syntax tree) */
typedef enum 
  {
    VAR,
    /* constants */
    CST_BOOL,CST_INT,
    /* data structures */
    LIST,DICT,
    /* operators */
    F_NEG,F_ABS,F_ADD,F_SUB,F_MUL,F_DIV,F_MOD,F_POW,
    F_IF,F_MIN,F_MAX,
    F_EQ,F_NE,F_GE,F_GT,F_LE,F_LT,
    F_NOT,F_AND,F_OR,F_XOR,F_IFF,
    /* symbolic constants */
    SYMB_NIL,SYMB_EQ,SYMB_NE,SYMB_GE,SYMB_GT,SYMB_LE,SYMB_LT,

    /* for internal use by the parser */
    _DICTKEY
  } NodeType;

#endif
