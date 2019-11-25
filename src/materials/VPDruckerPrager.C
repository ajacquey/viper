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

#include "VPDruckerPrager.h"

registerADMooseObject("ViperApp", VPDruckerPrager);

defineADValidParams(
    VPDruckerPrager,
    VPSingleVarUpdate,
    params.addClassDescription("Viscoplastic update based on a Drucker-Prager yield function.");
    params.addRequiredRangeCheckedParam<Real>("friction_angle",
                                              "friction_angle >= 0.0",
                                              "The friction angle for the Drucker-Prager yield.");
    params.addRangeCheckedParam<Real>("dilation_angle",
                                      0.0,
                                      "dilation_angle >= 0.0",
                                      "The dilation angle for the Drucker-Prager yield.");
    params.addRangeCheckedParam<Real>(
        "cohesion", 0.0, "cohesion >= 0.0", "The cohesion for the Drucker-Prager yield."););

template <ComputeStage compute_stage>
VPDruckerPrager<compute_stage>::VPDruckerPrager(const InputParameters & parameters)
  : VPSingleVarUpdate<compute_stage>(parameters),
    _phi(getParam<Real>("friction_angle")),
    _psi(getParam<Real>("dilation_angle")),
    _C(getParam<Real>("cohesion"))
{
  _alpha = std::sqrt(3.0) * std::sin(_phi * libMesh::pi / 180.0);
  _beta = std::sqrt(3.0) * std::sin(_psi * libMesh::pi / 180.0);
  _k = std::sqrt(3.0) * _C * std::cos(_phi * libMesh::pi / 180.0);
}

template <ComputeStage compute_stage>
ADReal
VPDruckerPrager<compute_stage>::yieldFunction(const ADReal & gamma_vp)
{
  return (_eqv_stress_tr - 3.0 * _G * gamma_vp * _dt) -
         _alpha * (_pressure_tr + _K * _beta * gamma_vp * _dt) - _k;
}

template <ComputeStage compute_stage>
ADReal
VPDruckerPrager<compute_stage>::yieldFunctionDeriv(const ADReal & /*gamma_vp*/)
{
  return -(3.0 * _G + _alpha * _beta * _K) * _dt;
}

template <ComputeStage compute_stage>
void
VPDruckerPrager<compute_stage>::preReturnMap()
{
  _pressure_tr = -_stress_tr.trace() / 3.0;
  _eqv_stress_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();
}

template <ComputeStage compute_stage>
void
VPDruckerPrager<compute_stage>::postReturnMap()
{
}

template <ComputeStage compute_stage>
ADRankTwoTensor
VPDruckerPrager<compute_stage>::reformPlasticStrainTensor(const ADReal & gamma_vp)
{
  ADRankTwoTensor flow_dir =
      (_eqv_stress_tr != 0.0) ? _stress_tr.deviatoric() / _eqv_stress_tr : ADRankTwoTensor();

  ADRankTwoTensor delta_gamma = 1.5 * gamma_vp * _dt * flow_dir;
  delta_gamma.addIa(_beta * gamma_vp * _dt / 3.0);

  return delta_gamma;
}
