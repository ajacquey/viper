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

#include "VPStressAuxBase.h"

class VPVonMisesStressAux;

template <>
InputParameters validParams<VPVonMisesStressAux>();

class VPVonMisesStressAux : public VPStressAuxBase
{
public:
  VPVonMisesStressAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
};