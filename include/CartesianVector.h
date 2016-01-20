  /// \file CaloHit.h
/*
 *
 * CaloHit.h header template automatically generated by a class generator
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

#ifndef CARTESIANVECTOR_H
#define CARTESIANVECTOR_H

#include <cmath>
#include <limits>
#include <ostream>

/**
 * @brief CartesianVector class
 */
class CartesianVector
{
public:

	/**
		*  @brief  Constructor, create a vector from the cartesian coordinates of the end point,
		*          origin at (0,0,0,)
		*
		*  @param  x the end point x coordinate
		*  @param  y the end point y coordinate
		*  @param  z the end point z coordinate
		*/
	CartesianVector(float x = 0.f, float y = 0.f, float z = 0.f);

	/**
		*  @brief  Copy constructor
		*
		*  @param  rhs the cartesian vector to copy
		*/
	CartesianVector(const CartesianVector &rhs);

	/**
		*  @brief  Set the values of cartesian vector components
		*
		*  @param  x the x coordinate
		*  @param  y the y coordinate
		*  @param  z the z coordinate
		*/
	void setValues(float x, float y, float z);

	/**
		*  @brief  Get the cartesian x coordinate
		*
		*  @return The cartesian x coordinate
		*/
	float getX() const;

	/**
		*  @brief  Get the cartesian y coordinate
		*
		*  @return The cartesian y coordinate
		*/
	float getY() const;

	/**
		*  @brief  Get the cartesian z coordinate
		*
		*  @return The cartesian z coordinate
		*/
	float getZ() const;

	/**
		*  @brief  Get the magnitude
		*
		*  @return The magnitude
		*/
	float getMagnitude() const;

	/**
		*  @brief  Get the magnitude squared
		*
		*  @return The magnitude squared
		*/
	float getMagnitudeSquared() const;

	/**
		*  @brief  Get the dot product of the cartesian vector with a second cartesian vector
		*
		*  @param  rhs the second cartesian vector
		*
		*  @return The dot product
		*/
	float getDotProduct(const CartesianVector &rhs) const;

	/**
		*  @brief  Get the cross product of the cartesian vector with a second cartesian vector
		*
		*  @param  rhs the second cartesian vector
		*
		*  @return The cross product
		*/
	CartesianVector getCrossProduct(const CartesianVector &rhs) const;

	/**
		*  @brief  Get the cosine of the opening angle of the cartesian vector with respect to a second cartesian vector
		*
		*  @param  rhs the second cartesian vector
		*
		*  @return The cosine of the opening angle
		*/
	float getCosOpeningAngle(const CartesianVector &rhs) const;

	/**
		*  @brief  Get the opening angle of the cartesian vector with respect to a second cartesian vector
		*
		*  @param  rhs the second cartesian vector
		*
		*  @return The opening angle
		*/
	float getOpeningAngle(const CartesianVector &rhs) const;

	/**
		*  @brief  Get the spherical coordinates of the cartesian vector
		*
		*  @param  radius the magnitude of the vector
		*  @param  phi the azimuth of the vector
		*  @param  theta the inclination of the vector
		*/
	void getSphericalCoordinates(float &radius, float &phi, float &theta) const;

	/**
		*  @brief  Get the cylindrical coordinates of the cartesian vector (x/y .. radius, z .. z)
		*
		*  @param  radius the radius (x,y-plane) of the vector
		*  @param  phi the azimuth of the vector
		*  @param  z the z position of the vector
		*/
	void getCylindricalCoordinates(float &radius, float &phi, float &z) const;

	/**
		*  @brief  Get a unit vector in the direction of the cartesian vector
		*
		*  @return The unit vector
		*/
	CartesianVector getUnitVector() const;

	/**
		*  @brief  Cartesian vector assignment operator
		*
		*  @param  rhs the cartesian vector to assign
		*/
	CartesianVector &operator=(const CartesianVector &rhs);

	/**
		*  @brief  Cartesian vector += operator
		*
		*  @param  rhs the cartesian vector to add
		*/
	CartesianVector &operator+=(const CartesianVector &rhs);

	/**
		*  @brief  Cartesian vector -= operator
		*
		*  @param  rhs the cartesian vector to subtract
		*/
	CartesianVector &operator-=(const CartesianVector &rhs);

	/**
		*  @brief  Cartesian vector *= operator
		*
		*  @param  scalar the scalar to multiply
		*/
	CartesianVector &operator*=(const double scalar);

	/**
		*  @brief  Cartesian vector == operator
		*
		*  @param  rhs the cartesian vector to compare
		*/
	bool operator==(const CartesianVector &rhs) const;

private:
    float   m_x;                ///< The x coordinate
    float   m_y;                ///< The y coordinate
    float   m_z;                ///< The z coordinate
};

/**
 *  @brief  Cartesian vector addition operator
 * 
 *  @param  lhs first cartesian vector, to which the second is added
 *  @param  rhs second cartesian vector, which is added to the first
 */
CartesianVector operator+(const CartesianVector &lhs, const CartesianVector &rhs);

/**
 *  @brief  Cartesian vector subtraction operator
 * 
 *  @param  lhs first cartesian vector, from which the second is subtracted
 *  @param  rhs second cartesian vector, which is subtracted from the first
 */
CartesianVector operator-(const CartesianVector &lhs, const CartesianVector &rhs);

/**
 *  @brief  Cartesian vector multiplication with scalar operator
 * 
 *  @param  lhs the cartesian vector to be multiplied by the scalar
 *  @param  scalar the value of the scalar
 */
CartesianVector operator*(const CartesianVector &lhs, const double scalar);

/**
 *  @brief  Operator to dump cartesian vector properties to an ostream
 *
 *  @param  stream the target ostream
 *  @param  cartesianVector the cartesian vector
 */
std::ostream &operator<<(std::ostream & stream, const CartesianVector& cartesianVector);

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector::CartesianVector(float x, float y, float z) :
    m_x(x),
    m_y(y),
    m_z(z)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector::CartesianVector(const CartesianVector &rhs) :
    m_x(rhs.m_x),
    m_y(rhs.m_y),
    m_z(rhs.m_z)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void CartesianVector::setValues(float x, float y, float z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getX() const
{
    return m_x;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getY() const
{
    return m_y;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getZ() const
{
    return m_z;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getMagnitude() const
{
    return std::sqrt(this->getMagnitudeSquared());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getMagnitudeSquared() const
{
    return ((m_x * m_x) + (m_y * m_y) + (m_z * m_z));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getDotProduct(const CartesianVector &rhs) const
{
    return ((m_x * rhs.m_x) + (m_y * rhs.m_y) + (m_z * rhs.m_z));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector CartesianVector::getCrossProduct(const CartesianVector &rhs) const
{
    return CartesianVector( (m_y * rhs.m_z) - (rhs.m_y * m_z),
                            (m_z * rhs.m_x) - (rhs.m_z * m_x),
                            (m_x * rhs.m_y) - (rhs.m_x * m_y));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float CartesianVector::getOpeningAngle(const CartesianVector &rhs) const
{
    return std::acos(this->getCosOpeningAngle(rhs));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector &CartesianVector::operator=(const CartesianVector &rhs)
{
    this->setValues(rhs.m_x, rhs.m_y, rhs.m_z);
    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector &CartesianVector::operator+=(const CartesianVector &rhs)
{
    this->setValues(m_x + rhs.m_x, m_y + rhs.m_y, m_z + rhs.m_z);
    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector &CartesianVector::operator-=(const CartesianVector &rhs)
{
    this->setValues(m_x - rhs.m_x, m_y - rhs.m_y, m_z - rhs.m_z);
    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector &CartesianVector::operator*=(const double scalar)
{
    this->setValues(static_cast<float>(m_x * scalar), static_cast<float>(m_y * scalar), static_cast<float>(m_z * scalar));
    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool CartesianVector::operator==(const CartesianVector &rhs) const
{
    return ( (std::fabs(m_x - rhs.m_x) < std::numeric_limits<float>::epsilon()) &&
        (std::fabs(m_y - rhs.m_y) < std::numeric_limits<float>::epsilon()) &&
        (std::fabs(m_z - rhs.m_z) < std::numeric_limits<float>::epsilon()) );
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector operator+(const CartesianVector &lhs, const CartesianVector &rhs)
{
    return CartesianVector(lhs.getX() + rhs.getX(), lhs.getY() + rhs.getY(), lhs.getZ() + rhs.getZ());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector operator-(const CartesianVector &lhs, const CartesianVector &rhs)
{
    return CartesianVector(lhs.getX() - rhs.getX(), lhs.getY() - rhs.getY(), lhs.getZ() - rhs.getZ());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CartesianVector operator*(const CartesianVector &lhs, const double scalar)
{
    return CartesianVector(static_cast<float>(lhs.getX() * scalar), static_cast<float>(lhs.getY() * scalar), static_cast<float>(lhs.getZ() * scalar));
}

#endif // #ifndef CARTESIANVECTOR_H
