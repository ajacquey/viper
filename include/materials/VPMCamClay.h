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

#include "VPTwoVarUpdate.h"

template <ComputeStage>
class VPMCamClay;

declareADValidParams(VPMCamClay);

template <ComputeStage compute_stage>
class VPMCamClay : public VPTwoVarUpdate<compute_stage>
{
public:
  VPMCamClay(const InputParameters & parameters);

protected:
  virtual ADReal yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d) override;
  virtual void overStress(const ADReal & gamma_v,
                          const ADReal & gamma_d,
                          ADReal & over_v,
                          ADReal & over_d) override;
  virtual void overStressDerivV(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_v,
                                ADReal & over_d_v) override;
  virtual void overStressDerivD(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_d,
                                ADReal & over_d_d) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap() override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_v,
                                                    const ADReal & gamma_d) override;
  virtual void
  calculateProjection(const ADReal & chi_v, const ADReal & chi_d, ADReal & chi_v0, ADReal & chi_d0);

  const Real _phi;
  const Real _pcr0;
  const Real _L;
  Real _alpha;

  ADReal _pressure_tr;
  ADReal _eqv_stress_tr;

  ADReal _pcr_tr;
  ADReal _A_tr;
  ADReal _B_tr;

  ADReal _chi_v_tr;
  ADReal _chi_d_tr;

  usingTwoVarUpdateMembers;
};