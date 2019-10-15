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

#include "VPViscoPlasticModel.h"

defineADValidParams(VPViscoPlasticModel,
                    ADMaterial,
                    params.addClassDescription("Base class for the viscoplastic correction.");
                    params.set<bool>("compute") = false;
                    params.suppressParameter<bool>("compute"););

template <ComputeStage compute_stage>
VPViscoPlasticModel<compute_stage>::VPViscoPlasticModel(const InputParameters & parameters)
  : ADMaterial<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
void
VPViscoPlasticModel<compute_stage>::setQp(unsigned int qp)
{
  _qp = qp;
}

adBaseClass(VPViscoPlasticModel);