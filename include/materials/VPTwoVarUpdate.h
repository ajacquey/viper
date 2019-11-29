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

#define usingTwoVarUpdateMembers                                                                   \
  usingViscoPlasticUpdateMembers;                                                                  \
  using VPTwoVarUpdate<compute_stage>::_stress_tr;                                                 \
  using VPTwoVarUpdate<compute_stage>::_K;                                                         \
  using VPTwoVarUpdate<compute_stage>::_G

template <ComputeStage>
class VPTwoVarUpdate;

declareADValidParams(VPTwoVarUpdate);

template <ComputeStage compute_stage>
class VPTwoVarUpdate : public VPViscoPlasticUpdate<compute_stage>
{
public:
  VPTwoVarUpdate(const InputParameters & parameters);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual void returnMap(ADReal & gamma_v, ADReal & gamma_d);
  virtual void
  residual(const ADReal & gamma_v, const ADReal & gamma_d, ADReal & resv, ADReal & resd);
  virtual void jacobian(const ADReal & gamma_v,
                        const ADReal & gamma_d,
                        ADReal & jacvv,
                        ADReal & jacdd,
                        ADReal & jacvd,
                        ADReal & jacdv);
  virtual ADReal yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d) = 0;
  virtual void
  overStress(const ADReal & gamma_v, const ADReal & gamma_d, ADReal & over_v, ADReal & over_d) = 0;
  virtual void overStressDerivV(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_v,
                                ADReal & over_d_v) = 0;
  virtual void overStressDerivD(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_d_v,
                                ADReal & over_d_d) = 0;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_v,
                                                    const ADReal & gamma_d) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d) = 0;

  ADRankTwoTensor _stress_tr;
  ADReal _K;
  ADReal _G;

  usingViscoPlasticUpdateMembers;
};