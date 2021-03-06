/* autogl_uims.c : user interface management system */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "autogl_math.h"
#include "autogl_uims.h"

#include "autogl_uims_inside.h"



/* Module Variables */

#define MAX_PROPERTYS 1000
static int NPropertys;
static AutoGL_Property Propertys[MAX_PROPERTYS];

static AutoGL_CallbackType ViewRedrawCallback;

static AutoGL_CallbackType ViewEventCallback;




/* Property */

static AutoGL_Property *GetNextProperty 
(AutoGL_PropertyType type, const char *name)
{
  AutoGL_Property *property;
  int nameLength;

  assert(NPropertys < MAX_PROPERTYS);
  property = &Propertys[NPropertys];
  NPropertys++;

  property->type = type;

  assert(name != NULL);
  nameLength = strlen (name);
  assert(nameLength < AUTOGL_PROPERTY_NAME_LENGTH);
  strcpy (property->name, name);

  strcpy (property->label, "");
  property->callback = NULL;
  property->booleanValue = NULL;
  property->stringValue = NULL;
  property->stringLength = 0;
  property->integerValue = NULL;
  property->minIntegerValue = 0;
  property->maxIntegerValue = 0;
  property->nItems = 0;
  property->realValue = NULL;
  property->minRealValue = 0.0;
  property->maxRealValue = 0.0;
  
  return property;
}

static AutoGL_Property *GetCurrentProperty (void)
{
  assert(0 < NPropertys);
  return &Propertys[NPropertys - 1];
}

void AutoGL_AddComment (void)
{
  GetNextProperty (AUTOGL_PROPERTY_COMMENT, "");
}

void AutoGL_AddCallback 
(AutoGL_CallbackType callback, 
 const char *name)
{
  AutoGL_Property *property;
  
  property = GetNextProperty (AUTOGL_PROPERTY_CALLBACK, name);
  assert(callback != NULL);
  property->callback = callback;
}

void AutoGL_AddBoolean 
(int *value_IO, 
 const char *name)
{
  AutoGL_Property *property;

  property = GetNextProperty (AUTOGL_PROPERTY_BOOLEAN, name);
  assert(value_IO != NULL);
  property->booleanValue = value_IO;
}

void AutoGL_AddString 
(char *str_IO, 
 const char *name, int length)
{
  AutoGL_Property *property;

  property = GetNextProperty (AUTOGL_PROPERTY_STRING, name);
  assert(str_IO != NULL);
  property->stringValue = str_IO;
  property->stringLength = length;
}

void AutoGL_AddInteger 
(int *value_IO, 
 const char *name)
{
  AutoGL_Property *property;

  property = GetNextProperty (AUTOGL_PROPERTY_INTEGER, name);
  assert(value_IO != NULL);
  property->integerValue = value_IO;
}

void AutoGL_SetIntegerRange 
(int minValue, int maxValue)
{
  AutoGL_Property *property;

  property = GetCurrentProperty ();
  assert(property->type == AUTOGL_PROPERTY_INTEGER);
  assert(minValue <= *property->integerValue);
  assert(*property->integerValue <= maxValue);
  property->minIntegerValue = minValue;
  property->maxIntegerValue = maxValue;
}

void AutoGL_AddIntegerItem 
(const char *label)
{
  AutoGL_Property *property;
  int labelLength;

  property = GetCurrentProperty ();
  assert(property->type == AUTOGL_PROPERTY_INTEGER);
  assert(property->nItems < AUTOGL_MAX_PROPERTY_ITEMS);
  assert(label != NULL);
  labelLength = strlen (label);
  assert(labelLength < AUTOGL_PROPERTY_LABEL_LENGTH);
  strcpy (property->itemLabels[property->nItems], label);
  property->nItems++;
}

void AutoGL_AddReal 
(double *value_IO, 
 const char *name)
{
  AutoGL_Property *property;

  property = GetNextProperty (AUTOGL_PROPERTY_REAL, name);
  assert(value_IO != NULL);
  property->realValue = value_IO;
}

void AutoGL_SetRealRange 
(double minValue, double maxValue)
{
  AutoGL_Property *property;

  property = GetCurrentProperty ();
  assert(property->type == AUTOGL_PROPERTY_REAL);
  assert(minValue <= *property->realValue);
  assert(*property->realValue <= maxValue);
  property->minRealValue = minValue;
  property->maxRealValue = maxValue;
}

void AutoGL_SetLabel 
(const char *label)
{
  AutoGL_Property *property;
  int labelLength;

  property = GetCurrentProperty ();
  labelLength = strlen (label);
  assert(labelLength < AUTOGL_PROPERTY_LABEL_LENGTH);
  strcpy (property->label, label);
}

int AutoGL_GetNPropertys (void)
{
  return NPropertys;
}

AutoGL_Property *AutoGL_GetProperty (int propertyId)
{
  assert(0 <= propertyId);
  assert(propertyId < NPropertys);
  return &Propertys[propertyId];
}

char *AutoGL_GetPropertyLabel 
(AutoGL_Property *property)
{
  char *label;

  label = property->name;
  if (0 < strlen (property->label)) {
    label = property->label;
  }

  return label;
}

void AutoGL_SetPropertyString 
(AutoGL_Property *property, 
 const char buf[])
{
  if (property->stringValue != NULL) {
    strcpy (property->stringValue, buf);
  } else if (property->integerValue != NULL) {
    sscanf (buf, "%d", property->integerValue);
  } else if (property->realValue != NULL) {
    sscanf (buf, "%lf", property->realValue);
  } else {
    assert (0);
  }
}

void AutoGL_GetPropertyString 
(AutoGL_Property *property, 
 char buf_OUT[256])
{
  if (property->stringValue != NULL) {
    strcpy (buf_OUT, property->stringValue);
  } else if (property->integerValue != NULL) {
    sprintf (buf_OUT, "%d", *property->integerValue);
  } else if (property->realValue != NULL) {
    sprintf (buf_OUT, "%f", *property->realValue);
  } else {
    assert (0);
  }
}

void AutoGL_CallCallback 
(AutoGL_Property *property)
{
  AutoGL_CallbackType callback = property->callback;
  
  (*callback) ();
}



/* View Callback */

void AutoGL_SetViewRedrawCallback 
(AutoGL_CallbackType callback)
{
  assert(callback != NULL);
  ViewRedrawCallback = callback;
}

void AutoGL_SetViewEventCallback 
(AutoGL_CallbackType callback)
{
  assert(callback != NULL || callback == NULL);
  ViewEventCallback = callback;
}

AutoGL_CallbackType AutoGL_GetViewRedrawCallback (void)
{
  return ViewRedrawCallback;
}

AutoGL_CallbackType AutoGL_GetViewEventCallback (void)
{
  return ViewEventCallback;
}

void AutoGL_CallViewRedrawCallback (void)
{
  if (ViewRedrawCallback != NULL) {
    (*ViewRedrawCallback) ();
  }
}

void AutoGL_CallViewEventCallback (void)
{
  if (ViewEventCallback != NULL) {
    (*ViewEventCallback) ();
  }
}





