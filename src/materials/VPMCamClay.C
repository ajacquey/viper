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

#include "VPMCamClay.h"

registerADMooseObject("ViperApp", VPMCamClay);

defineADValidParams(
    VPMCamClay,
    VPTwoVarUpdate,
    params.addClassDescription("Viscoplastic update based on a modified Cam-Clay yield function.");
    params.addRequiredRangeCheckedParam<Real>("friction_angle",
                                              "friction_angle > 0.0",
                                              "The friction angle for the critical state line.");
    params.addRequiredRangeCheckedParam<Real>(
        "critical_pressure",
        "critical_pressure > 0.0",
        "The critical pressure of the modified Cam-Clay yield.");
    params.addRangeCheckedParam<Real>("critical_pressure_hardening",
                                      0.0,
                                      "critical_pressure_hardening >= 0.0",
                                      "The hardening parameter in the exponential for the critical "
                                      "pressure of the modified Cam-Clay yield."););

template <ComputeStage compute_stage>
VPMCamClay<compute_stage>::VPMCamClay(const InputParameters & parameters)
  : VPTwoVarUpdate<compute_stage>(parameters),
    _phi(getParam<Real>("friction_angle")),
    _pcr0(getParam<Real>("critical_pressure")),
    _L(getParam<Real>("critical_pressure_hardening"))
{
  _alpha = std::sqrt(3.0) * std::sin(_phi * libMesh::pi / 180.0);
}

template <ComputeStage compute_stage>
ADReal
VPMCamClay<compute_stage>::yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d)
{
  // Here we calculate the yield function in the dissipative stress space:
  // chi_v = pressure - 0.5 * pc
  // chi_d = eqv_stress
  ADReal chi_v =
      _chi_v_tr - _K * gamma_v * _dt; // + 0.5 * _pc_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  ADReal chi_d = _chi_d_tr - 3.0 * _G * gamma_d * _dt;
  ADReal A = _A_tr;
  ADReal B = _B_tr;
  return Utility::pow<2>(chi_v / A) + Utility::pow<2>(chi_d / B) - 1.0;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::overStress(const ADReal & gamma_v,
                                      const ADReal & gamma_d,
                                      ADReal & over_v,
                                      ADReal & over_d)
{
  ADReal chi_v = _chi_v_tr - _K * gamma_v * _dt;
  // + 0.5 * _pc_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  ADReal chi_d = _chi_d_tr - 3.0 * _G * gamma_d * _dt;

  ADReal chi_v0 = 0.0, chi_d0 = 0.0;
  calculateProjection(chi_v, chi_d, chi_v0, chi_d0);

  over_v = chi_v - chi_v0;
  over_d = chi_d - chi_d0;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::overStressDerivV(const ADReal & gamma_v,
                                            const ADReal & gamma_d,
                                            ADReal & over_v_v,
                                            ADReal & over_d_v)
{
  // How to calculate the derivatives of the projections?
  ADReal dchi_v = -_K * _dt;
  // + 0.5 * _pc_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  ADReal dchi_d = 0.0;

  over_v_v = dchi_v;
  over_d_v = dchi_d;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::overStressDerivD(const ADReal & gamma_v,
                                            const ADReal & gamma_d,
                                            ADReal & over_v_d,
                                            ADReal & over_d_d)
{
  // How to calculate the derivatives of the projections?
  ADReal dchi_v = 0.0;
  // + 0.5 * _pc_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  ADReal dchi_d = -3.0 * _G * _dt;

  over_v_d = dchi_v;
  over_d_d = dchi_d;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::preReturnMap()
{
  _pressure_tr = -_stress_tr.trace() / 3.0;
  _eqv_stress_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();

  _pcr_tr = _pcr0; // * std::exp(-_L * _plastic_strain_old[_qp].trace());
  _A_tr = 0.5 * _pcr_tr;
  _B_tr = _alpha * 0.5 * _pcr_tr;

  _chi_v_tr = _pressure_tr - 0.5 * _pcr_tr;
  _chi_d_tr = _eqv_stress_tr;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::postReturnMap()
{
}

template <ComputeStage compute_stage>
ADRankTwoTensor
VPMCamClay<compute_stage>::reformPlasticStrainTensor(const ADReal & gamma_v, const ADReal & gamma_d)
{
  ADRankTwoTensor flow_dir =
      (_eqv_stress_tr != 0.0) ? _stress_tr.deviatoric() / _eqv_stress_tr : ADRankTwoTensor();

  ADRankTwoTensor delta_gamma = 1.5 * gamma_d * _dt * flow_dir;
  delta_gamma.addIa(-gamma_v * _dt / 3.0);

  return delta_gamma;
}

template <ComputeStage compute_stage>
void
VPMCamClay<compute_stage>::calculateProjection(const ADReal & chi_v,
                                               const ADReal & chi_d,
                                               ADReal & chi_v0,
                                               ADReal & chi_d0)
{
  // Directions
  ADReal rho_tr = std::sqrt(Utility::pow<2>(chi_v) + Utility::pow<2>(chi_d));
  ADReal ev = (rho_tr != 0.0) ? chi_v / rho_tr : 0.0;
  ADReal ed = (rho_tr != 0.0) ? chi_d / rho_tr : 0.0;
  ADReal A = _A_tr;
  ADReal B = _B_tr;

  ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / A) + Utility::pow<2>(ed / B)));

  chi_v0 = rho_0 * ev;
  chi_d0 = rho_0 * ed;
}