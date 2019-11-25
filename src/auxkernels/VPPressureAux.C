/******************************************************************************/
/*                            This file is part of                            */
/*                       VIPER, a MOOSE-based application                     */
/*                          VIsco-Plastic bEnchmaRks                          */
/*                                                                            */
/*  Copyright (C) 2019 by Antoine B. Jacquey, Thomas Poulet and Mustafa Sari  */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*   CSIRO, The Commonwealth Scientific and Industrial Research Organisation  */
/*                                                                            */
/*            Licensed under GNU Lesser General Public License v2.1           */
/*                       please see LICENSE for details                       */
/*                 or http://www.gnu.org/licenses/lgpl.html                   */
/******************************************************************************/

#include "VPPressureAux.h"

registerMooseObject("ViperApp", VPPressureAux);

template <>
InputParameters
validParams<VPPressureAux>()
{
  InputParameters params = validParams<VPStressAuxBase>();
  params.addClassDescription("Calculates the pressure.");
  return params;
}

VPPressureAux::VPPressureAux(const InputParameters & parameters) : VPStressAuxBase(parameters) {}

Real
VPPressureAux::computeValue()
{
  return -_stress[_qp].trace() / 3.0;
}