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

#pragma once

#include "ADMaterial.h"

template <ComputeStage>
class VPViscoPlasticModel;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(VPViscoPlasticModel);

template <ComputeStage compute_stage>
class VPViscoPlasticModel : public ADMaterial<compute_stage>
{
public:
  VPViscoPlasticModel(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const ADReal & K,
                                  const ADReal & G,
                                  ADRankTwoTensor & elastic_strain_incr) = 0;
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  usingMaterialMembers;
};