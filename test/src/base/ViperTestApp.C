//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ViperTestApp.h"
#include "ViperApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

template <>
InputParameters
validParams<ViperTestApp>()
{
  InputParameters params = validParams<ViperApp>();
  return params;
}

ViperTestApp::ViperTestApp(InputParameters parameters) : MooseApp(parameters)
{
  ViperTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ViperTestApp::~ViperTestApp() {}

void
ViperTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ViperApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ViperTestApp"});
    Registry::registerActionsTo(af, {"ViperTestApp"});
  }
}

void
ViperTestApp::registerApps()
{
  registerApp(ViperApp);
  registerApp(ViperTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ViperTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ViperTestApp::registerAll(f, af, s);
}
extern "C" void
ViperTestApp__registerApps()
{
  ViperTestApp::registerApps();
}
