/***************************************************************************
                          qgsnumeric.h  -  description
                             -------------------
    begin                : October 2012
    copyright            : (C) 2012 by Hugo Mercier
    email                : hugo dot mercier at oslandia dot com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef QGSNUMERIC_H
#define QGSNUMERIC_H

#include <limits>

/**
 * shortcut to get NaN for double
 */
inline double NaN() { return std::numeric_limits<double>::quiet_NaN(); }

/**
 * shortcut to test NaN for double
 */
inline bool isNaN( const double& value ) { return value != value; }

#endif
