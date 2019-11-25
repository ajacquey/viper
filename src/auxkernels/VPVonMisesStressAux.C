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

#include "VPVonMisesStressAux.h"

registerMooseObject("ViperApp", VPVonMisesStressAux);

template <>
InputParameters
validParams<VPVonMisesStressAux>()
{
  InputParameters params = validParams<VPStressAuxBase>();
  params.addClassDescription("Calculates the Von Mises stress.");
  return params;
}

VPVonMisesStressAux::VPVonMisesStressAux(const InputParameters & parameters)
  : VPStressAuxBase(parameters)
{
}

Real
VPVonMisesStressAux::computeValue()
{
  RankTwoTensor stress_dev = _stress[_qp].deviatoric();
  return std::sqrt(1.5) * stress_dev.L2norm();
}
