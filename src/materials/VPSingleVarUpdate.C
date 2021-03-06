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

#include "VPSingleVarUpdate.h"
#include "ElasticityTensorTools.h"

defineADValidParams(
    VPSingleVarUpdate,
    VPViscoPlasticUpdate,
    params.addClassDescription("Base class for a single variable viscoplastic update."););

template <ComputeStage compute_stage>
VPSingleVarUpdate<compute_stage>::VPSingleVarUpdate(const InputParameters & parameters)
  : VPViscoPlasticUpdate<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
void
VPSingleVarUpdate<compute_stage>::viscoPlasticUpdate(ADRankTwoTensor & stress,
                                                     const RankFourTensor & Cijkl,
                                                     ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with a single variable (usually scalar viscoplastic strain rate)
  // We are trying to find the zero of the function F which is defined as:
  // F(gamma_vp) = yield - eta * gamma_vp^(1 / n)
  // gamma_vp: scalar viscoplastic strain rate (scalar)
  // yield: the yield function
  // eta: the viscoplastic viscosity
  // n: exponent for Perzyna-like flow rule
  // flow rule: gamma_vp = (yield / eta)^n

  // Trial stress
  _stress_tr = stress;
  // Elastic moduli
  _K = ElasticityTensorTools::getIsotropicBulkModulus(Cijkl);
  _G = ElasticityTensorTools::getIsotropicShearModulus(Cijkl);

  // Initialize plastic strain increment
  _plastic_strain_incr[_qp].zero();

  // Pre return map calculations (model specific)
  preReturnMap();

  // Check yield function
  _yield_function[_qp] = yieldFunction(0.0);
  if (_yield_function[_qp] <= _abs_tol) // Elastic
    return;

  // Viscoplastic update
  ADReal gamma_vp = returnMap();

  // Update quantities
  _yield_function[_qp] = yieldFunction(gamma_vp);
  _plastic_strain_incr[_qp] = reformPlasticStrainTensor(gamma_vp);
  elastic_strain_incr -= _plastic_strain_incr[_qp];
  stress -= Cijkl * _plastic_strain_incr[_qp];
  postReturnMap();
}

template <ComputeStage compute_stage>
ADReal
VPSingleVarUpdate<compute_stage>::returnMap()
{
  // Initialize scalar viscoplastic strain rate
  ADReal gamma_vp = 0.0;

  // Initial residual
  ADReal res_ini = residual(gamma_vp);

  ADReal res = res_ini;
  ADReal jac = jacobian(gamma_vp);

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    gamma_vp -= res / jac;

    res = residual(gamma_vp);
    jac = jacobian(gamma_vp);

    // Convergence check
    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
      return gamma_vp;
  }
  throw MooseException("VPSingleVarUpdate: maximum number of iterations exceeded in 'returnMap'!");
}

template <ComputeStage compute_stage>
ADReal
VPSingleVarUpdate<compute_stage>::residual(const ADReal & gamma_vp)
{
  ADReal res = yieldFunction(gamma_vp);
  if (gamma_vp != 0.0)
    res -= _eta_p * std::pow(gamma_vp, 1.0 / _n);

  return res;
}

template <ComputeStage compute_stage>
ADReal
VPSingleVarUpdate<compute_stage>::jacobian(const ADReal & gamma_vp)
{
  ADReal jac = yieldFunctionDeriv(gamma_vp);
  if (gamma_vp != 0.0)
    jac -= _eta_p / _n * std::pow(gamma_vp, 1.0 / _n - 1.0);

  return jac;
}

adBaseClass(VPSingleVarUpdate);