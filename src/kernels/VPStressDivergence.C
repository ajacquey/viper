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

#include "VPStressDivergence.h"
#include "libmesh/quadrature.h"

registerADMooseObject("ViperApp", VPStressDivergence);

defineADValidParams(
    VPStressDivergence, ADKernel, params.addClassDescription("Solid momentum kernel.");
    params.set<bool>("use_displaced_mesh") = false;
    params.addRequiredParam<unsigned int>("component",
                                          "An integer corresponding to the direction "
                                          "the variable this kernel acts in (0 for x, "
                                          "1 for y, 2 for z)."););

template <ComputeStage compute_stage>
VPStressDivergence<compute_stage>::VPStressDivergence(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _component(getParam<unsigned int>("component")),
    _stress(getADMaterialProperty<RankTwoTensor>("stress"))
{
}

template <ComputeStage compute_stage>
ADReal
VPStressDivergence<compute_stage>::computeQpResidual()
{
  return _stress[_qp].row(_component) * _grad_test[_i][_qp];
}

adBaseClass(VPStressDivergence);