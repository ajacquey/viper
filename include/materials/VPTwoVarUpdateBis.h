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

#define usingTwoVarUpdateBisMembers                                                                   \
  usingViscoPlasticUpdateMembers;                                                                  \
  using VPTwoVarUpdateBis<compute_stage>::_stress_tr;                                                 \
  using VPTwoVarUpdateBis<compute_stage>::_K;                                                         \
  using VPTwoVarUpdateBis<compute_stage>::_G

template <ComputeStage>
class VPTwoVarUpdateBis;

declareADValidParams(VPTwoVarUpdateBis);

template <ComputeStage compute_stage>
class VPTwoVarUpdateBis : public VPViscoPlasticUpdate<compute_stage>
{
public:
  VPTwoVarUpdateBis(const InputParameters & parameters);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual void returnMap(ADReal & p, ADReal & q);
  virtual void residual(const ADReal & p,
                        const ADReal & q,
                        const ADReal & p_y,
                        const ADReal & q_y,
                        ADReal & resv,
                        ADReal & resd) = 0;
  virtual void jacobian(const ADReal & p,
                        const ADReal & q,
                        const ADReal & p_y,
                        const ADReal & q_y,
                        ADReal & jacvv,
                        ADReal & jacdd,
                        ADReal & jacvd,
                        ADReal & jacdv) = 0;
  virtual ADReal yieldFunction(const ADReal & p, const ADReal & q) = 0;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & p,
                                                    const ADReal & q) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap(const ADReal & p, const ADReal & q) = 0;
  virtual void calculateProjection(const ADReal & p,
                                   const ADReal & q,
                                   ADReal & p_y,
                                   ADReal & q_y) = 0;

  ADRankTwoTensor _stress_tr;
  ADReal _K;
  ADReal _G;

  usingViscoPlasticUpdateMembers;
};
