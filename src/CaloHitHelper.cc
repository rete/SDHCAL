  /// \file CaloHitHelper.cc
/*
 *
 * CaloHitHelper.cc source template automatically generated by a class generator
 * Creation date : sam. mai 31 2014
 *
 * This file is part of SDHCAL libraries.
 * 
 * SDHCAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * SDHCAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SDHCAL.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Remi Ete
 * @copyright CNRS , IPNL
 */


#include "CaloHitHelper.h"

// sdhcal
#include "CaloHit.h"
#include "CartesianVector.h"

// lcio
#include "EVENT/CalorimeterHit.h"
#include "UTIL/CellIDDecoder.h"

// std
#include <algorithm>

void CaloHitHelper::loadCaloHitCollection(const EVENT::LCCollection *pCollection,
		CaloHitList *pCaloHitList,
		OrderedCaloHitList *pOrderedCaloHitList)
{
	UTIL::CellIDDecoder<EVENT::CalorimeterHit> decoder(pCollection);

	const std::string codingString(pCollection->getParameters().getStringVal(LCIO::CellIDEncoding));
	const std::string iDecoderStr("I");
	const std::string jDecoderStr("J");
	const std::string kDecoderStr((codingString.find("K-1") == std::string::npos) ? "K" : "K-1");

	for(unsigned int elt=0 ; elt<pCollection->getNumberOfElements() ; elt++)
	{
		EVENT::CalorimeterHit *pCalorimeterHit = dynamic_cast<EVENT::CalorimeterHit *>(pCollection->getElementAt(elt));

		CartesianVector position(
				pCalorimeterHit->getPosition()[0],
				pCalorimeterHit->getPosition()[1],
				pCalorimeterHit->getPosition()[2]);

		CaloHitCell cell;
		cell.m_iCell = decoder(pCalorimeterHit)[iDecoderStr.c_str()];
		cell.m_jCell = decoder(pCalorimeterHit)[jDecoderStr.c_str()];
		cell.m_layer = decoder(pCalorimeterHit)[kDecoderStr.c_str()];

		SemiDigitalThreshold threshold;
		const float energy = pCalorimeterHit->getEnergy();
		if(energy == 1.f)
		{
			threshold = THRESHOLD_1;
		}
		else if(energy == 2.f)
		{
			threshold = THRESHOLD_2;
		}
		else if(energy == 3.f)
		{
			threshold = THRESHOLD_3;
		}
		else
		{
			std::cerr << "Couldn't load calorimeter hit collection. \n"
					"Bad energy conversion for SDHCAL threshold. Expected is 1.f, 2.f or 3.f. Given was " << energy << std::endl;
			return;
		}

		// build the calo hit and fill the collection(s)
		CaloHit *pCaloHit = new CaloHit(position, cell, threshold);
		pCaloHitList->insert(pCaloHit);

		if(NULL != pOrderedCaloHitList)
		{
		 (*pOrderedCaloHitList)[cell.m_layer].insert(pCaloHit);
		}
	}
}

//------------------------------------------------------------------------------------------------

void CaloHitHelper::calculateCaloHitDensities(const OrderedCaloHitList &orderedCaloHitList)
{
	for(OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin() , layerEndIter = orderedCaloHitList.end() ;
			layerEndIter != layerIter ; ++layerIter)
	{
		Layer layer = layerIter->first;

		for(CaloHitList::iterator iter = layerIter->second.begin() , endIter = layerIter->second.end() ; endIter != iter ; ++iter)
		{
			CaloHit *pCaloHit = *iter;

			CaloHitHelper::calculateCaloHitDensity2D(pCaloHit, layerIter->second);
		}
	}

	for(OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin() , layerEndIter = orderedCaloHitList.end() ;
			layerEndIter != layerIter ; ++layerIter)
	{
		Layer layer = layerIter->first;

		for(CaloHitList::iterator iter = layerIter->second.begin() , endIter = layerIter->second.end() ; endIter != iter ; ++iter)
		{
			CaloHit *pCaloHit = *iter;

			CaloHitHelper::calculateCaloHitDensity3D(pCaloHit, orderedCaloHitList);
		}
	}
}

//------------------------------------------------------------------------------------------------

void CaloHitHelper::calculateCaloHitDensity2D(CaloHit *pCaloHit, const CaloHitList &caloHitList)
{
	for(CaloHitList::iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
	{
		CaloHit *pOtherCaloHit = *iter;

		if(abs(static_cast<int>(pOtherCaloHit->getCell().m_iCell) - static_cast<int>(pCaloHit->getCell().m_iCell)) < 2
		&& abs(static_cast<int>(pOtherCaloHit->getCell().m_jCell) - static_cast<int>(pCaloHit->getCell().m_jCell)) < 2
		&& pOtherCaloHit->getCell().m_layer == pCaloHit->getCell().m_layer)
		{
			pCaloHit->getDensityInfo().m_weightedDensity2D += pOtherCaloHit->getSemiDigitalThreshold();
			pCaloHit->getDensityInfo().m_unWeightedDensity2D ++;
		}
	}
}

//------------------------------------------------------------------------------------------------

void CaloHitHelper::calculateCaloHitDensity3D(CaloHit *pCaloHit, const OrderedCaloHitList &orderedCaloHitList)
{
	const Layer &layer = pCaloHit->getCell().m_layer;
	OrderedCaloHitList::const_iterator previousLayerIter = orderedCaloHitList.end();
	OrderedCaloHitList::const_iterator nextLayerIter = orderedCaloHitList.end();

	if(layer != 0)
	{
		previousLayerIter = orderedCaloHitList.find(layer-1);
	}

	nextLayerIter = orderedCaloHitList.find(layer+1);

	if(orderedCaloHitList.end() != previousLayerIter)
	{
		for(CaloHitList::iterator iter = previousLayerIter->second.begin() , endIter = previousLayerIter->second.end() ;
				endIter != iter ; ++iter)
		{
			CaloHit *pOtherCaloHit = *iter;

			if(abs(static_cast<int>(pOtherCaloHit->getCell().m_iCell) - static_cast<int>(pCaloHit->getCell().m_iCell)) < 2
			&& abs(static_cast<int>(pOtherCaloHit->getCell().m_jCell) - static_cast<int>(pCaloHit->getCell().m_jCell)) < 2 )
			{
				pCaloHit->getDensityInfo().m_weightedDensity3D += pOtherCaloHit->getSemiDigitalThreshold();
				pCaloHit->getDensityInfo().m_unWeightedDensity3D++;
			}
		}
	}

	if(orderedCaloHitList.end() != nextLayerIter)
	{
		for(CaloHitList::iterator iter = nextLayerIter->second.begin() , endIter = nextLayerIter->second.end() ;
				endIter != iter ; ++iter)
		{
			CaloHit *pOtherCaloHit = *iter;

			if(abs(static_cast<int>(pOtherCaloHit->getCell().m_iCell) - static_cast<int>(pCaloHit->getCell().m_iCell)) < 2
			&& abs(static_cast<int>(pOtherCaloHit->getCell().m_jCell) - static_cast<int>(pCaloHit->getCell().m_jCell)) < 2 )
			{
				pCaloHit->getDensityInfo().m_weightedDensity3D += pOtherCaloHit->getSemiDigitalThreshold();
				pCaloHit->getDensityInfo().m_unWeightedDensity3D ++;
			}
		}
	}

	pCaloHit->getDensityInfo().m_weightedDensity3D += pCaloHit->getDensityInfo().m_weightedDensity2D;
	pCaloHit->getDensityInfo().m_unWeightedDensity3D += pCaloHit->getDensityInfo().m_unWeightedDensity2D;
}



void CaloHitHelper::getNThreshold1(const CaloHitList &caloHitList, unsigned int &nHit1)
{
	nHit1 = 0;

	for(CaloHitList::iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
	{
		CaloHit *pCaloHit = *iter;

		if(THRESHOLD_1 == pCaloHit->getSemiDigitalThreshold())
			nHit1 ++;
	}
}

void CaloHitHelper::getNThreshold2(const CaloHitList &caloHitList, unsigned int &nHit2)
{
	nHit2 = 0;

	for(CaloHitList::iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
	{
		CaloHit *pCaloHit = *iter;

		if(THRESHOLD_2 == pCaloHit->getSemiDigitalThreshold())
			nHit2 ++;
	}
}

void CaloHitHelper::getNThreshold3(const CaloHitList &caloHitList, unsigned int &nHit3)
{
	nHit3 = 0;

	for(CaloHitList::iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
	{
		CaloHit *pCaloHit = *iter;

		if(THRESHOLD_3 == pCaloHit->getSemiDigitalThreshold())
			nHit3 ++;
	}
}

void CaloHitHelper::getNThresholds(const CaloHitList &caloHitList,
		unsigned int &nHit1, unsigned int &nHit2, unsigned int &nHit3)
{
	nHit1 = 0;
	nHit2 = 0;
	nHit3 = 0;

	for(CaloHitList::iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
	{
		CaloHit *pCaloHit = *iter;

		if(THRESHOLD_1 == pCaloHit->getSemiDigitalThreshold())
			nHit1 ++;
		if(THRESHOLD_2 == pCaloHit->getSemiDigitalThreshold())
			nHit2 ++;
		if(THRESHOLD_3 == pCaloHit->getSemiDigitalThreshold())
			nHit3 ++;
	}
}

