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

#include "VPVolStrainRateAux.h"

registerMooseObject("ViperApp", VPVolStrainRateAux);

template <>
InputParameters
validParams<VPVolStrainRateAux>()
{
  InputParameters params = validParams<VPStrainAuxBase>();
  params.addClassDescription("Calculates the volumetric strain rate of the given tensor.");
  return params;
}

VPVolStrainRateAux::VPVolStrainRateAux(const InputParameters & parameters)
  : VPStrainAuxBase(parameters)
{
}

Real
VPVolStrainRateAux::computeValue()
{
  return (*_strain_incr)[_qp].trace() / _dt;
}
