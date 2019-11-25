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

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class VPStrainAuxBase;

template <>
InputParameters validParams<VPStrainAuxBase>();

class VPStrainAuxBase : public DerivativeMaterialInterface<AuxKernel>
{
public:
  VPStrainAuxBase(const InputParameters & parameters);
  static MooseEnum strainType();

protected:
  MooseEnum _strain_type;
  std::string _strain_name;
  const MaterialProperty<RankTwoTensor> * _strain_incr;
};