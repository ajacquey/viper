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

#include "ADMaterial.h"

template <ComputeStage>
class VPMechMaterial;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;
template <ComputeStage>
class VPViscoPlasticUpdate;

declareADValidParams(VPMechMaterial);

template <ComputeStage compute_stage>
class VPMechMaterial : public ADMaterial<compute_stage>
{
public:
  VPMechMaterial(const InputParameters & parameters);
  void initialSetup() override;
  void displacementIntegrityCheck();

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  virtual void computeQpStrainIncrement();
  virtual void computeQpSmallStrain(const ADRankTwoTensor & grad_tensor,
                                    const RankTwoTensor & grad_tensor_old);
  virtual void computeQpFiniteStrain(const ADRankTwoTensor & grad_tensor,
                                     const RankTwoTensor & grad_tensor_old);
  virtual void computeQpElasticityTensor();
  virtual void computeQpStress();
  virtual ADRankTwoTensor spinRotation(const ADRankTwoTensor & tensor);

  // Coupled variables
  const unsigned int _ndisp;
  std::vector<const ADVariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_old;

  // Strain parameters
  const unsigned int _strain_model;

  // Elastic parameters
  const Real _bulk_modulus;
  const Real _shear_modulus;

  // Viscoplastoc model
  const bool _has_vp;

  // Strain properties
  ADMaterialProperty(RankTwoTensor) & _strain_increment;
  ADMaterialProperty(RankTwoTensor) & _spin_increment;
  ADMaterialProperty(RankTwoTensor) & _elastic_strain_incr;

  // Stress properties
  ADMaterialProperty(RankTwoTensor) & _stress;
  const MaterialProperty<RankTwoTensor> & _stress_old;

  // Viscoplastic model
  VPViscoPlasticUpdate<compute_stage> * _vp_model;

  // Elasticity tensor
  RankFourTensor _Cijkl;

  usingMaterialMembers;
};