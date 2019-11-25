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

#include "VPEqvStrainRateAux.h"

registerMooseObject("ViperApp", VPEqvStrainRateAux);

template <>
InputParameters
validParams<VPEqvStrainRateAux>()
{
  InputParameters params = validParams<VPStrainAuxBase>();
  params.addClassDescription("Calculates the equivalent strain rate of the given tensor.");
  return params;
}

VPEqvStrainRateAux::VPEqvStrainRateAux(const InputParameters & parameters)
  : VPStrainAuxBase(parameters)
{
}

Real
VPEqvStrainRateAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * (*_strain_incr)[_qp].deviatoric().L2norm() / _dt;
}
