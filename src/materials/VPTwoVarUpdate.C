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

#include "VPTwoVarUpdate.h"
#include "ElasticityTensorTools.h"

defineADValidParams(
    VPTwoVarUpdate,
    VPViscoPlasticUpdate,
    params.addClassDescription("Base class for a single variable viscoplastic update."););

template <ComputeStage compute_stage>
VPTwoVarUpdate<compute_stage>::VPTwoVarUpdate(const InputParameters & parameters)
  : VPViscoPlasticUpdate<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdate<compute_stage>::viscoPlasticUpdate(ADRankTwoTensor & stress,
                                                  const RankFourTensor & Cijkl,
                                                  ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with two variables (usually scalar volumetric and deviatoric
  // viscoplastic strain rates)
  // We are trying to find the zero of the two residual functions Rv, Rd which are defined as:
  // Rv(gamma_v, gamma_d) = p - p0 - eta * gamma_v^(1 / n)
  // Rd(gamma_v, gamma_d) = q - q0 - eta * gamma_d^(1 / n)
  // gamma_v: scalar volumetric viscoplastic strain rate
  // gamma_d: scalar deviatoric viscoplastic strain rate
  // p0, q0: the projection of the pressure and deviatoric stress on the yield function
  // eta: the viscoplastic viscosity
  // n: exponent for Perzyna-like flow rule

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
  _yield_function[_qp] = yieldFunction(0.0, 0.0);
  if (_yield_function[_qp] <= _abs_tol) // Elastic
    return;

  // Viscoplastic update
  ADReal gamma_v = 0.0, gamma_d = 0.0;
  returnMap(gamma_v, gamma_d);

  // Update quantities
  _yield_function[_qp] = yieldFunction(gamma_v, gamma_d);
  _plastic_strain_incr[_qp] = reformPlasticStrainTensor(gamma_v, gamma_d);
  elastic_strain_incr -= _plastic_strain_incr[_qp];
  stress -= Cijkl * _plastic_strain_incr[_qp];
  postReturnMap(gamma_v, gamma_d);
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdate<compute_stage>::returnMap(ADReal & gamma_v, ADReal & gamma_d)
{
  // Initial residual
  ADReal resv_ini = 0.0, resd_ini = 0.0;
  residual(0.0, 0.0, resv_ini, resd_ini);
  ADReal res_ini = std::sqrt(Utility::pow<2>(resv_ini) + Utility::pow<2>(resd_ini));
  ADReal resv = resv_ini, resd = resd_ini;
  ADReal res = res_ini;

  // Initial jacobian
  ADReal jacvv = 0.0, jacdd = 0.0, jacvd = 0.0, jacdv = 0.0;
  jacobian(0.0, 0.0, jacvv, jacdd, jacvd, jacdv);

  // Useful stuff
  ADReal jac_full = jacvv * jacdd - jacvd * jacdv;
  ADReal resv_full = jacdd * resv - jacvd * resd;
  ADReal resd_full = jacvv * resd - jacdv * resv;

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    gamma_v -= resv_full / jac_full;
    gamma_d -= resd_full / jac_full;

    residual(gamma_v, gamma_d, resv, resd);
    jacobian(gamma_v, gamma_d, jacvv, jacdd, jacvd, jacdv);
    jac_full = jacvv * jacdd - jacvd * jacdv;
    resv_full = jacdd * resv - jacvd * resd;
    resd_full = jacvv * resd - jacdv * resv;
    res = std::sqrt(Utility::pow<2>(resv) + Utility::pow<2>(resd));

    // Convergence check
    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
      return;
  }
  throw MooseException(
      "VPTwoVarUpdate: maximum number of iterations exceeded in 'returnMap'!\n Residuals: ",
      res,
      "\n Vol strain rate: ",
      gamma_v,
      "\n Dev strain Rate: ",
      gamma_d,
      "\n");
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdate<compute_stage>::residual(const ADReal & gamma_v,
                                        const ADReal & gamma_d,
                                        ADReal & resv,
                                        ADReal & resd)
{
  overStress(gamma_v, gamma_d, resv, resd);
  if (gamma_v != 0.0)
    resv -= _eta_p * std::pow(gamma_v, 1.0 / _n);
  if (gamma_d != 0.0)
    resd -= _eta_p * std::pow(gamma_d, 1.0 / _n);
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdate<compute_stage>::jacobian(const ADReal & gamma_v,
                                        const ADReal & gamma_d,
                                        ADReal & jacvv,
                                        ADReal & jacdd,
                                        ADReal & jacvd,
                                        ADReal & jacdv)
{
  overStressDerivV(gamma_v, gamma_d, jacvv, jacdv);
  overStressDerivD(gamma_v, gamma_d, jacvd, jacdd);
  if (gamma_v != 0.0)
    jacvv -= _eta_p / _n * std::pow(gamma_v, 1.0 / _n - 1.0);
  if (gamma_d != 0.0)
    jacdd -= _eta_p / _n * std::pow(gamma_d, 1.0 / _n - 1.0);
}

adBaseClass(VPTwoVarUpdate);