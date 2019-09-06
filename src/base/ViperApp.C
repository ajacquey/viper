#include "ViperApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

template <>
InputParameters
validParams<ViperApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

ViperApp::ViperApp(InputParameters parameters) : MooseApp(parameters)
{
  ViperApp::registerAll(_factory, _action_factory, _syntax);
}

ViperApp::~ViperApp() {}

void
ViperApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"ViperApp"});
  Registry::registerActionsTo(af, {"ViperApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ViperApp::registerApps()
{
  registerApp(ViperApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ViperApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ViperApp::registerAll(f, af, s);
}
extern "C" void
ViperApp__registerApps()
{
  ViperApp::registerApps();
}
