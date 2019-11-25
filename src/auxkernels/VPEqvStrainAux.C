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

#include "VPEqvStrainAux.h"

registerMooseObject("ViperApp", VPEqvStrainAux);

template <>
InputParameters
validParams<VPEqvStrainAux>()
{
  InputParameters params = validParams<VPStrainAuxBase>();
  params.addClassDescription("Calculates the equivalent strain of the given tensor.");
  return params;
}

VPEqvStrainAux::VPEqvStrainAux(const InputParameters & parameters) : VPStrainAuxBase(parameters) {}

Real
VPEqvStrainAux::computeValue()
{
  return _u_old[_qp] + std::sqrt(2.0 / 3.0) * (*_strain_incr)[_qp].deviatoric().L2norm();
}
