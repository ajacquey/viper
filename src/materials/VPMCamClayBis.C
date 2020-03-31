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
        "The critical pressure of the modified Cam-Clay yield.");
    MooseEnum projection_method("ray_intersection=1 gradient_potential=2", "ray_intersection");
    params.addParam<MooseEnum>("projection_method", projection_method,
        "The type of stress projection to the yield enveloppe.");
    params.addRangeCheckedParam<unsigned int>(
        "projection_max_iterations",
        10,
        "projection_max_iterations >= 1",
        "The maximum number of iterations for the iterative projection to "\
        "the yield enveloppe");
    params.addRangeCheckedParam<Real>("projection_abs_tolerance",
        1.0e-14,
        "projection_abs_tolerance > 0.0",
        "The absolute tolerance for the iterative projection."););

template <ComputeStage compute_stage>
VPMCamClayBis<compute_stage>::VPMCamClayBis(const InputParameters & parameters)
  : VPTwoVarUpdateBis<compute_stage>(parameters),
    _M(getParam<Real>("critical_state_line_slope")),
    _pc(-getParam<Real>("critical_pressure")),
    _eps_dot_0(std::pow(_eta_p, -_n)),
    _projection_method(getParam<MooseEnum>("projection_method")),
    _proj_abs_tol(getParam<Real>("projection_abs_tolerance")),
    _proj_max_its(getParam<unsigned int>("projection_max_iterations"))
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
                                       const ADReal & p_y,
                                       const ADReal & q_y,
                                       ADReal & resv,
                                       ADReal & resd)
{
  ADReal gamma_v, gamma_d;
  gamma_v = std::abs(_p_tr - p) / (_K * _dt);
  gamma_d = std::abs(_q_tr - q) / (3 * _G * _dt);

  resv = std::abs(p - p_y);
  if (gamma_v != 0.0)
    resv -= _eta_p * std::pow(gamma_v, 1.0 / _n);
  resd = std::abs(q - q_y);
  if (gamma_d != 0.0)
    resd -= _eta_p * std::pow(gamma_d, 1.0 / _n);
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::jacobian(const ADReal & p,
                                       const ADReal & q,
                                       const ADReal & p_y,
                                       const ADReal & q_y,
                                       ADReal & jacvv,
                                       ADReal & jacdd,
                                       ADReal & jacvd,
                                       ADReal & jacdv)
{
  // First, compute yield point derivatives
  ADReal dpy_dp = 0.0;
  ADReal dpy_dq = 0.0;
  ADReal dqy_dp = 0.0;
  ADReal dqy_dq = 0.0;
  switch (_projection_method)
  {
    case 1: // Ray intersection from centre of ellipse
    {
      ADReal x5 = (p - 0.5 * _pc);
      ADReal x1 = std::abs(_pc) / (2 * std::pow(Utility::pow<2>(q / _M) + x5 * x5, 1.5));
      dpy_dp = x1 * q * q / (_M * _M);
      dpy_dq = -x1 * q * x5 / (_M * _M);
      dqy_dp = x1 * q * x5;
      dqy_dq = x1 * x5 * x5;
      break;
    }

    case 2: // Follow gradient of potential
    {
      if (q > 0)
      {
        if (p == 0.5 * _pc)
        {
          dpy_dp = std::pow(2 * _M * std::abs(_pc) / q, Utility::pow<2>(_M));
        } else {
          // General case,
          ADReal M2 = Utility::pow<2>(_M);
          ADReal pc2 = Utility::pow<2>(_pc);
          ADReal x1 = p - 0.5 * _pc;
          ADReal x2 = std::exp(-2 * _t_y);
          ADReal x3 = std::exp(-2 * _t_y / M2);
          ADReal x4 = Utility::pow<2>(2 * x1 * x2);
          ADReal x5 = x1 * Utility::pow<2>(2 * x2) / (pc2 + x4 * (M2 - 1.0));
          ADReal x6 = (pc2 - x4) / (pc2  + x4 *  (M2 - 1.0));
          dpy_dp = x2 * (1 - M2 * x1 * x5);
          dqy_dq = x3 * (1 - x6);
          dpy_dq = -M2 * x1 * x2 * x6 / q;
          dpy_dq = -q * x3 * x5;
        }
      } else {
        // q == 0
        dqy_dq = std::pow(std::abs(_pc / (2 * p - _pc)), 1.0 / Utility::pow<2>(_M));
      }
      break;
    }

    default:
      mooseError("Unknown projection type. Specify 'ray_intersection' "\
                 "or 'potential_gradient'!");
  }

  // Second, compute jacobians
  ADReal x1 = 3 * _G * _eps_dot_0 * _dt;
  ADReal x2 = std::abs(_q_tr - q) / x1;
  if (x2 > 0.0)
    x2 = std::pow(x2, 1.0 / _n - 1.0);
  jacdd = (1.0 - dqy_dq) * (q > q_y ? 1 : -1) + x2 * (_q_tr > q ? 1 : -1) / (_n * x1);
  jacdv = -dqy_dp * (q > q_y ? 1 : -1);

  ADReal x3 = _K * _eps_dot_0 * _dt;
  ADReal x4 = std::abs(_p_tr - p) / x3;
  if (x4 > 0.0)
    x4 = std::pow(x4, 1.0 / _n - 1.0);
  jacvv = (1.0 - dpy_dp) * (p > p_y ? 1 : -1) + x4  * (_p_tr > p ? 1 : -1) / (_n * x3);
  jacvd = -dpy_dq * (p > p_y ? 1 : -1);
}

template <ComputeStage compute_stage>
void
VPMCamClayBis<compute_stage>::preReturnMap()
{
  _p_tr = _stress_tr.trace() / 3.0;
  _q_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();
  _t_y = -1.0;
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
  switch (_projection_method)
  {
    case 1: // Ray intersection from centre of ellipse
      q_y = 0.5 * std::abs(_pc) * q /
        std::sqrt(Utility::pow<2>(q / _M) + Utility::pow<2>(p - 0.5 * _pc));
      if (q > 0)
        p_y = 0.5 * _pc + (p - 0.5 * _pc) * q_y / q;
      else if (p < 0.5 *_pc)
        p_y = _pc;
      else
        p_y = 0.0;
      break;

    case 2: // Follow gradient of potential
      // TODO, case M==1 => same as ray_intersection to save time (?)
      if (q > 0)
      {
        if (p == 0.5 * _pc)
        {
          p_y = 0.5 * _pc;
          q_y = _M * std::abs(p_y);
        } else {
          // Generic case, using Newton Raphson
          //std::cout << "Starting NR for p=" << p << ", q=" << q << ", pc=" << _pc << ", M=" << _M << std::endl;
          ADReal M2 = Utility::pow<2>(_M);
          ADReal x1 = Utility::pow<2>(q / _M);
          ADReal x2 = Utility::pow<2>(p - 0.5 * _pc);
          ADReal x3 = Utility::pow<2>(0.5 * _pc);
          ADReal x4 = -4 * Utility::pow<2>(q / M2);
          _t_y = 0.0;
          ADReal phi = x1 * std::exp(-4 * _t_y / M2) + x2 * std::exp(-4 * _t_y) - x3;
          ADReal phi_prime = x4 * std::exp(-4 * _t_y / M2) - 4 * x2 * std::exp(-4 * _t_y);
          unsigned int iter = 0;
          while ((std::fabs(phi) > _proj_abs_tol) && (iter < _proj_max_its))
          {
              iter += 1;
              //std::cout << "iter=" << iter << ", phi=" << phi << std::endl;
              _t_y -= phi / phi_prime;
              phi = x1 * std::exp(-4 * _t_y / M2) + x2 * std::exp(-4 * _t_y) - x3;
              phi_prime = x4 * std::exp(-4 * _t_y / M2) - 4 * x2 * std::exp(-4 * _t_y);
          }
          if (iter == _proj_max_its)
          {
            throw MooseException(
              "VPMCamClayBis: maximum number of iterations exceeded in "\
              "'calculateProjection'!\n iter:", iter, ", phi: ", phi, "\n");
          }
          p_y = 0.5 * _pc + (p - 0.5 * _pc) * std::exp(-2 * _t_y);
          q_y = q * std::exp(-2 * _t_y / M2);
          //std::cout << "Newton Raphson converged in " << iter << " iterations, t_y="
          //    << _t_y << ", phi=" << phi << ", p_y=" << p_y << ", q_y=" << q_y <<std::endl;
        }
      } else if (p < 0.5 *_pc) {
        // q == 0, p <= _pc
        p_y = _pc;
        q_y = 0.0;
      } else {
        // q == 0, p >= 0
        p_y = 0.0;
        q_y = 0.0;
      }
      break;

    default:
      mooseError("Unknown projection type. Specify 'ray_intersection' "\
                 "or 'potential_gradient'!");
  }
}

