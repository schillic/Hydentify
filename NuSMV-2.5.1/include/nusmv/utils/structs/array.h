/**CHeaderFile*****************************************************************

  FileName    [array.h]

  PackageName [utils.structs]

  Synopsis    [Macro for declaring structures based on dynamic arrays]

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

  Revision    [$Id: array.h,v 1.1.2.1 2008-09-15 07:58:47 nusmv Exp $]

******************************************************************************/


#ifndef __STRUCT_ARRAY_H__
#define __STRUCT_ARRAY_H__

#include "utils/array.h"

#define DECLARE_ARRAY_H(name, eltype)    \
                                         \
typedef void (* name##_DESTROY)(eltype); \
array_t* name##_access();                \
eltype name##_get(int pos);              \
void name##_set(int pos, eltype val);    \
void name##_append(eltype val);          \
int name##_get_length()

#define DECLARE_ARRAY_C(name, eltype)        \
                                             \
static void name##_init();                   \
static void name##_quit();                   \
static void name##_quit_fun(name##_DESTROY); \
array_t* name##_access();                    \
eltype name##_get(int pos);                  \
void name##_set(int pos, eltype val);        \
void name##_append(eltype val);              \
int name##_get_length()



#define INSTANTIATE_ARRAY(name, eltype) \
                                                \
static array_t* name = NULL;                    \
static void name##_init()                       \
{                                               \
  nusmv_assert(name == NULL);                   \
  name = array_alloc(eltype, 1);                \
  nusmv_assert(name != NULL);                   \
}                                               \
static void name##_quit()                       \
{                                               \
  nusmv_assert(name != NULL);                   \
  array_free(name);                             \
  name = NULL;                                  \
}                                               \
static void name##_quit_fun(name##_DESTROY f)   \
{                                               \
  int i; eltype el;                             \
  nusmv_assert(name != NULL);                   \
  arrayForEachItem(eltype,name,i,el) {          \
    if (el != NULL) f(el);                      \
  }                                             \
  array_free(name);                             \
  name = NULL;                                  \
}                                               \
array_t* name##_access()                        \
{                                               \
  nusmv_assert(name != NULL);                   \
  return name;                                  \
}                                               \
eltype name##_get(int pos)                      \
{                                               \
  nusmv_assert(name != NULL);                   \
  return array_fetch(eltype, name, pos);        \
}                                               \
void name##_set(int pos, eltype val)            \
{                                               \
  nusmv_assert(name != NULL);                   \
  array_insert(eltype, name, pos, val);         \
}                                               \
void name##_append(eltype val)                  \
{                                               \
  nusmv_assert(name != NULL);                   \
  array_insert_last(eltype, name, val);         \
}                                               \
int name##_get_length()                         \
{                                               \
  nusmv_assert(name != NULL);                   \
  return array_n(name);                         \
}                                               \


#endif /* __STRUCT_ARRAY_H__ */
