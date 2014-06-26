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
#ifndef _C_AST_h_
#define _C_AST_h_

#include "XMLParser_constants.h"

/**
 * 
 * a simple representation of an Abstract Syntax Tree (AST) for C
 *
 */

/*
 * abtract node (the field type is present in every node of the
 * tree). All other nodes in the tree "inherit" from this type.
 *
 * This is also the representation of a symbolic constant.
 */
typedef
struct 
{
  NodeType type;
} C_AST;


/*
 * representation of a variable (necessarily a leaf in the tree)
 */
typedef
struct 
{
  NodeType type;
  char *varName;
  int idVar;
} C_AST_VarNode;

/*
 * representation of a constant (boolean or integer -- necessarily a
 * leaf in the tree)
 */
typedef
struct 
{
  NodeType type;
  int val;
} C_AST_CstNode;

/*
 * representation of a function 
 */
typedef
struct 
{
  NodeType type;
  int nbarg;
  C_AST *args[];
} C_AST_FxNode;

/*
 * representation of a list 
 */
typedef
struct 
{
  NodeType type;
  int size; /* number of items in the list */
  C_AST *items[];
} C_AST_ListNode;

/*
 * entry in a dictionary
 */
typedef
struct 
{
  char *key;
  C_AST *value;
} C_AST_DictEntry;

/*
 * representation of a dictionary 
 */
typedef
struct 
{
  NodeType type;
  int size; /* number of items in the dictionary */
  C_AST_DictEntry items[];
} C_AST_DictNode;



#endif
