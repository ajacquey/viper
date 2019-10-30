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

#include "VPViscoPlasticUpdate.h"

#define usingSingleVarUpdateMembers                                                                \
  usingViscoPlasticUpdateMembers;                                                                  \
  using VPSingleVarUpdate<compute_stage>::_stress_tr;                                              \
  using VPSingleVarUpdate<compute_stage>::_K;                                                      \
  using VPSingleVarUpdate<compute_stage>::_G

template <ComputeStage>
class VPSingleVarUpdate;

declareADValidParams(VPSingleVarUpdate);

template <ComputeStage compute_stage>
class VPSingleVarUpdate : public VPViscoPlasticUpdate<compute_stage>
{
public:
  VPSingleVarUpdate(const InputParameters & parameters);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual ADReal returnMap();
  virtual ADReal residual(const ADReal & gamma_vp);
  virtual ADReal jacobian(const ADReal & gamma_vp);
  virtual ADReal yieldFunction(const ADReal & gamma_vp) = 0;
  virtual ADReal yieldFunctionDeriv(const ADReal & gamma_vp) = 0;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_vp) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap() = 0;

  ADRankTwoTensor _stress_tr;
  ADReal _K;
  ADReal _G;

  usingViscoPlasticUpdateMembers;
};