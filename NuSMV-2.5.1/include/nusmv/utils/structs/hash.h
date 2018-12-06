/**CHeaderFile*****************************************************************

  FileName    [array.h]

  PackageName [utils.structs]

  Synopsis    [Macro for declaring structures based on hash tables]

  Description []

  SeeAlso     []

  Author      [Roberto Cavada]

  Copyright   [
  This file is part of the ``utils.structs'' package of NuSMV version 2. 
  Copyright (C) 2008 by FBK-irst. 

  NuSMV version 2 is free software; you can redistribute it and/or 
  modify it under the terms of the GNU Lesser General Public 
  License as published by the Free Software Foundation; either 
  version 2 of the License, or (at your option) any later version.

  NuSMV version 2 is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public 
  License along with this library; if not, write to the Free Software 
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

  For more information on NuSMV see <http://nusmv.fbk.eu>
  or email to <nusmv-users@fbk.eu>.
  Please report bugs to <nusmv-users@fbk.eu>.

  To contact the NuSMV development board, email to <nusmv@fbk.eu>. ]

  Revision    [$Id: hash.h,v 1.1.2.2 2008-09-18 11:19:19 nusmv Exp $]

******************************************************************************/

#ifndef __STRUCTS_HASH_H__
#define __STRUCTS_HASH_H__


#include "node/node.h"
#include "utils/assoc.h"


/**
 * Declare and instantiate an hash table whose name is given, and
 * provides some functions to accesss the hash: 
 *  - name_init and name_quit 
 *  - name_lookup and name_insert
 */
/* -----------------------------------------*/
#define DECLARE_HASH_H(name, keytype, eltype) \
/* -----------------------------------------*/\
typedef void (* name##_DESTROY)(eltype);      \
eltype name##_lookup(keytype key);            \
void name##_insert(keytype key, eltype val)


/* -----------------------------------------*/
#define DECLARE_HASH_C(name, keytype, eltype) \
/* -----------------------------------------*/\
static void name##_init();                    \
static void name##_quit();                    \
static void name##_quit_fun(name##_DESTROY);  \
eltype name##_lookup(keytype key);            \
void name##_insert(keytype key, eltype val)


/* -----------------------------------------*/
#define DECLARE_HASH_HC(name, keytype, eltype) \
/* -----------------------------------------*/ \
typedef void (* name##_DESTROY)(eltype);       \
static void name##_init();                     \
static void name##_quit();                     \
static void name##_quit_fun(name##_DESTROY);   \
eltype name##_lookup(keytype key);             \
void name##_insert(keytype key, eltype val)


/* --------------------------------------------*/
#define INSTANTIATE_HASH(name, keytype, eltype)                      \
/* --------------------------------------------*/                    \
static hash_ptr name = NULL;                                         \
static void name##_init()                                            \
{                                                                    \
  nusmv_assert(name == NULL);                                        \
  name = new_assoc();                                                \
  nusmv_assert(name != NULL);                                        \
}                                                                    \
static void name##_quit()                                            \
{                                                                    \
  nusmv_assert(name != NULL);                                        \
  free_assoc(name);                                                  \
  name = NULL;                                                       \
}                                                                    \
static enum st_retval name##_quit_fun_aux(char* k, char* e, char* a) \
{                                                                    \
  name##_DESTROY f = (name##_DESTROY) a;                             \
  f((eltype) e);                                                     \
  return ASSOC_DELETE;                                               \
}                                                                    \
static void name##_quit_fun(name##_DESTROY fun)                      \
{                                                                    \
  nusmv_assert(name != NULL);                                        \
  clear_assoc_and_free_entries_arg(name, name##_quit_fun_aux,        \
                                   (char*) fun);                     \
  free_assoc(name);                                                  \
  name = NULL;                                                       \
}                                                                    \
eltype name##_lookup(keytype key)                                    \
{                                                                    \
  nusmv_assert(name != NULL);                                        \
  return (eltype) find_assoc(name, (node_ptr) key);                  \
}                                                                    \
void name##_insert(keytype key, eltype val)                          \
{                                                                    \
  nusmv_assert(name != NULL);                                        \
  insert_assoc(name, (node_ptr) key, (node_ptr) val);                \
}                                                                    \


#endif /* __STRUCTS_HASH_H__ */
