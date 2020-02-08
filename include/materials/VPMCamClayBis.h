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

#include "VPTwoVarUpdateBis.h"

template <ComputeStage>
class VPMCamClayBis;

declareADValidParams(VPMCamClayBis);

template <ComputeStage compute_stage>
class VPMCamClayBis : public VPTwoVarUpdateBis<compute_stage>
{
public:
  VPMCamClayBis(const InputParameters & parameters);

protected:
  virtual ADReal yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d) override;
  virtual void residual(const ADReal & p,
                        const ADReal & q,
                        ADReal & resv,
                        ADReal & resd) override;
  virtual void jacobian(const ADReal & p,
                        const ADReal & q,
                        ADReal & jacvv,
                        ADReal & jacdd,
                        ADReal & jacvd,
                        ADReal & jacdv) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap(const ADReal & p, const ADReal & q) override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & p,
                                                    const ADReal & q) override;
  virtual void
  calculateProjection(const ADReal & p, const ADReal & q, ADReal & p_y, ADReal & q_y);

  const Real _M;
  const Real _pc;
  const Real _eps_dot_0;
  const unsigned int _projection_method;
  const Real _proj_abs_tol; // for gradient_potential method only
  const unsigned int _proj_max_its; // for gradient_potential method only

  ADReal _p_tr;
  ADReal _q_tr;
  ADReal _t_y;

  usingTwoVarUpdateBisMembers;
};
