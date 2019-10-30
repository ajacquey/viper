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

#include "VPSingleVarUpdate.h"

template <ComputeStage>
class VPDruckerPrager;

declareADValidParams(VPDruckerPrager);

template <ComputeStage compute_stage>
class VPDruckerPrager : public VPSingleVarUpdate<compute_stage>
{
public:
  VPDruckerPrager(const InputParameters & parameters);

protected:
  virtual ADReal yieldFunction(const ADReal & gamma_vp) override;
  virtual ADReal yieldFunctionDeriv(const ADReal & gamma_vp) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap() override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_vp) override;

  Real _phi;
  Real _psi;
  Real _C;
  Real _alpha;
  Real _beta;
  Real _k;

  ADReal _pressure_tr;
  ADReal _eqv_stress_tr;

  usingSingleVarUpdateMembers;
};