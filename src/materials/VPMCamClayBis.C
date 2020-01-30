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

#include "VPMCamClayBis.h"

registerADMooseObject("ViperApp", VPMCamClayBis);

defineADValidParams(
    VPMCamClayBis,
    VPTwoVarUpdateBis,
    params.addClassDescription("Viscoplastic update based on a modified Cam-Clay yield function (bis).");
    params.addRequiredRangeCheckedParam<Real>("critical_state_line_slope",
                                              "critical_state_line_slope > 0.0",
                                              "The slope of the critical state line.");
    params.addRequiredRangeCheckedParam<Real>(
        "critical_pressure",
        "critical_pressure > 0.0",
        "The critical pressure of the modified Cam-Clay yield."););

template <ComputeStage compute_stage>
VPMCamClayBis<compute_stage>::VPMCamClayBis(const InputParameters & parameters)
  : VPTwoVarUpdateBis<compute_stage>(parameters),
    _M(getParam<Real>("critical_state_line_slope")),
    _pc(-getParam<Real>("critical_pressure")),
    _eps_dot_0(std::pow(_eta_p, -_n))
{
}

template <ComputeStage compute_stage>
ADReal
VPMCamClayBis<compute_stage>::yieldFunction(const ADReal & p, const ADReal & q)
{
  return Utility::pow<2>(q / _M) + p * (p - _pc);
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::residual(const ADReal & p,
                                       const ADReal & q,
                                       ADReal & resv,
                                       ADReal & resd)
{
  ADReal p_y = p, q_y = q;
  calculateProjection(p, q, p_y, q_y);

  ADReal gamma_v, gamma_d;
  if (q > 0)
  {
    gamma_v =  _M * _M * (_q_tr / q - 1.0) * (2 * p - _pc) / (6 * _G * _dt);
    gamma_d = (_q_tr - q) / (3 * _G * _dt);
  } else {
    // q ==  0
    gamma_v = (_p_tr - p) / (_K * _dt);
    gamma_d = 0;
  }
  resv = p - p_y - _eta_p * std::pow(gamma_v, 1.0 / _n);
  resd = q - q_y - _eta_p * std::pow(gamma_d, 1.0 / _n);
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::jacobian(const ADReal & p,
                                       const ADReal & q,
                                       ADReal & jacvv,
                                       ADReal & jacdd,
                                       ADReal & jacvd,
                                       ADReal & jacdv)
{
  if (q > 0)
  {
    ADReal x5 = (p - 0.5 * _pc);
    ADReal x1 = std::abs(_pc) / (2 * std::pow(Utility::pow<2>(q / _M) + x5 * x5, 1.5));
    ADReal x6 = 3 * _G * _eps_dot_0 * _dt;
    ADReal x2 = _M * _M / (_n * x6);
    ADReal x3 = _q_tr / q - 1.0;
    ADReal x4 = _n * x2 * x5 * x3;
    if ((x4 == 0.0) || (_n == 1)) // somehow std::pow(0, 0) returns (1,{nan,nan,nan,...
      x4 = 1.0;
    else
      x4 = std::pow(x4, 1.0 / _n - 1.0);
    ADReal x7 = (_q_tr - q) / x6;
    if ((x7 == 0.0) || (_n == 1))
      x7 = 1.0;
    else
      x7 = std::pow(x7, 1.0 / _n - 1.0);

    jacvv = 1.0 - x1 * q * q / (_M * _M) - x2 * x3 * x4;
    jacdd = 1.0 - x1 * x5 * x5 + x7 / (_n * x6);
    jacvd = x1 * q * x5 / (_M * _M) + x2 * x5 * _q_tr * x4 / (q * q);
    jacdv = -x1 * q * x5;
  } else {
    // q == 0
    jacvv = 1.0;
    jacdd = 1.0;
    jacvd = 0.0;
    jacdv = 0.0;
  }
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::preReturnMap()
{
  _p_tr = _stress_tr.trace() / 3.0;
  _q_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::postReturnMap(const ADReal & /*p*/, const ADReal & /*q*/)
{
}

template <ComputeStage compute_stage>
ADRankTwoTensor
VPMCamClayBis<compute_stage>::reformPlasticStrainTensor(const ADReal & p, const ADReal & q)
{
  ADRankTwoTensor plastic_strain_incr = ADRankTwoTensor();
  if (q > 0)
  {
    plastic_strain_incr.addIa(Utility::pow<2>(_M) * (2 * p - _pc) / 9.0);
    plastic_strain_incr += _stress_tr.deviatoric() * q / _q_tr;
    plastic_strain_incr *= (_q_tr / q -  1.0) / (2.0 * _G);
  } else {
    // q == 0
    plastic_strain_incr.addIa((p - _p_tr) / (3.0 * _K));
  }
  return plastic_strain_incr;
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::calculateProjection(const ADReal & p,
                                                  const ADReal & q,
                                                  ADReal & p_y,
                                                  ADReal & q_y)
{
  q_y = 0.5 * std::abs(_pc) * q /
    std::sqrt(Utility::pow<2>(q / _M) + Utility::pow<2>(p - 0.5 * _pc));
  if (q > 0)
    p_y = 0.5 * _pc + (p - 0.5 * _pc) * q_y / q;
  else if (p < 0.5 *_pc)
    p_y = _pc;
  else
    p_y = 0.0;
}
