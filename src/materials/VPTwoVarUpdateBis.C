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

#include "VPTwoVarUpdateBis.h"
#include "ElasticityTensorTools.h"

defineADValidParams(
    VPTwoVarUpdateBis,
    VPViscoPlasticUpdate,
    params.addClassDescription("Base class for a single variable viscoplastic update."););

template <ComputeStage compute_stage>
VPTwoVarUpdateBis<compute_stage>::VPTwoVarUpdateBis(const InputParameters & parameters)
  : VPViscoPlasticUpdate<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdateBis<compute_stage>::viscoPlasticUpdate(ADRankTwoTensor & stress,
                                                    const RankFourTensor & Cijkl,
                                                    ADRankTwoTensor & elastic_strain_incr)
{
  // Here we do an iterative update with two variables (usually p and q)
  // We are trying to find the zero of the two residual functions Rv, Rd which are defined as:
  // Rv(p, q) = p - p0 - eta * gamma_v^(1 / n)
  // Rd(p, q) = q - q0 - eta * gamma_d^(1 / n)
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
  ADReal p = _stress_tr.trace() / 3.0;
  ADReal q = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();
  _yield_function[_qp] = yieldFunction(p, q);
  if (_yield_function[_qp] <= _abs_tol) // Elastic
    return;

  // Viscoplastic update
  returnMap(p, q);

  // Update quantities
  _yield_function[_qp] = yieldFunction(p, q);
  _plastic_strain_incr[_qp] = reformPlasticStrainTensor(p, q);

  elastic_strain_incr -= _plastic_strain_incr[_qp];
  stress -= Cijkl * _plastic_strain_incr[_qp];
  postReturnMap(p, q);
}

template <ComputeStage compute_stage>
void
VPTwoVarUpdateBis<compute_stage>::returnMap(ADReal & p, ADReal & q)
{
  //std::cout << "\tEnter RMap, p="  << p << ",q="  << q << std::endl;

  // Initial residual
  ADReal resv_ini = 0.0, resd_ini = 0.0;
  residual(p, q, resv_ini, resd_ini);
  ADReal res_ini = std::sqrt(Utility::pow<2>(resv_ini) + Utility::pow<2>(resd_ini));
  ADReal resv = resv_ini, resd = resd_ini;
  ADReal res = res_ini;

  // Initial jacobian
  ADReal jacvv = 0.0, jacdd = 0.0, jacvd = 0.0, jacdv = 0.0;
  jacobian(p, q, jacvv, jacdd, jacvd, jacdv);

  // Useful stuff
  ADReal jac_full = jacvv * jacdd - jacvd * jacdv;
  ADReal resv_full = jacdd * resv - jacvd * resd;
  ADReal resd_full = jacvv * resd - jacdv * resv;

  // Newton loop
  for (unsigned int iter = 0; iter < _max_its; ++iter)
  {
    p -= resv_full / jac_full;
    q -= resd_full / jac_full;
    residual(p, q, resv, resd);
    jacobian(p, q, jacvv, jacdd, jacvd, jacdv);
    jac_full = jacvv * jacdd - jacvd * jacdv;
    resv_full = jacdd * resv - jacvd * resd;
    resd_full = jacvv * resd - jacdv * resv;
    res = std::sqrt(Utility::pow<2>(resv) + Utility::pow<2>(resd));

    /*std::cout << "\tNR, iter=" << iter ;
    ADReal abs_tol = std::abs(res);
    ADReal rel_tol = std::abs(res / res_ini);
    std::cout << ", p="  << p << ",q="  << q << ", abs_tol=" << abs_tol <<
        ", rel_tol=" << rel_tol << std::endl;*/

    if ((std::abs(res) <= _abs_tol) || (std::abs(res / res_ini) <= _rel_tol))
    {
      //std::cout << "\tNR CV, iter=" << iter << std::endl;
      return;
    }
  }
  throw MooseException(
      "VPTwoVarUpdateBis: maximum number of iterations exceeded in 'returnMap'!\n Residuals: ",
      res,
      "\n p: ",
      p,
      "\n q: ",
      q,
      "\n");
}

adBaseClass(VPTwoVarUpdateBis);
