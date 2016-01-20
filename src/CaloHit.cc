  /// \file CaloHit.cc
/*
 *
 * CaloHit.cc source template automatically generated by a class generator
 * Creation date : lun. avr. 28 2014
 *
 * This file is part of SDHCALEventDisplay libraries.
 * 
 * SDHCALEventDisplay is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * SDHCALEventDisplay is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SDHCALEventDisplay.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Remi Ete
 * @copyright CNRS , IPNL
 */


#include "CaloHit.h"

CaloHit::CaloHit(CartesianVector position, CaloHitCell cell, SemiDigitalThreshold semiDigitalThreshold) :
   m_position(position),
   m_cell(cell),
   m_semiDigitalThreshold(semiDigitalThreshold)
{

}

//-------------------------------------------------------------------------------------------

CaloHit::CaloHit(CaloHit *pCaloHit)
{
	m_position             = pCaloHit->m_position;
	m_cell                 = pCaloHit->m_cell;
	m_semiDigitalThreshold = pCaloHit->m_semiDigitalThreshold;
	m_densityInfo.m_weightedDensity2D = pCaloHit->m_densityInfo.m_weightedDensity2D;
	m_densityInfo.m_unWeightedDensity2D = pCaloHit->m_densityInfo.m_unWeightedDensity2D;
	m_densityInfo.m_weightedDensity3D = pCaloHit->m_densityInfo.m_weightedDensity3D;
	m_densityInfo.m_unWeightedDensity3D = pCaloHit->m_densityInfo.m_unWeightedDensity3D;
}

//-------------------------------------------------------------------------------------------

CaloHit::~CaloHit() 
{

}
