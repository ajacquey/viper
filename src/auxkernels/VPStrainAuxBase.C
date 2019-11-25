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

#include "VPStrainAuxBase.h"

template <>
InputParameters
validParams<VPStrainAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Base class for outputting strain values.");
  params.addParam<MooseEnum>("strain_type",
                             VPStrainAuxBase::strainType() = "total",
                             "The type of the strain tensor to output.");
  return params;
}

VPStrainAuxBase::VPStrainAuxBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _strain_type(getParam<MooseEnum>("strain_type"))
{
  switch (_strain_type)
  {
    case 1:
      _strain_name = "strain_increment";
      break;
    case 2:
      _strain_name = "elastic_strain_increment";
      break;
    case 3:
      _strain_name = "plastic_strain_increment";
      break;
    default:
      mooseError("VPStrainAuxBase: unknown strain type!");
  }
  _strain_incr = &getDefaultMaterialProperty<RankTwoTensor>(_strain_name);
}

MooseEnum
VPStrainAuxBase::strainType()
{
  return MooseEnum("total=1 elastic=2 plastic=3");
}
