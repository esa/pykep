/*****************************************************************************
 *   Copyright (C) 2004-2012 The PyKEP development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://keptoolbox.sourceforge.net/index.html                            *
 *   http://keptoolbox.sourceforge.net/credits.html                          *
 *                                                                           *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#ifndef ARRAY3D_OPERATIONS_H
#define ARRAY3D_OPERATIONS_H

#include <cmath>
#include"../astro_constants.h"

namespace kep_toolbox {
    template<class vettore3D>
    inline void sum(vettore3D& out, const vettore3D& v1,const vettore3D& v2){
        out[0] = v1[0]+v2[0];
        out[1] = v1[1]+v2[1];
        out[2] = v1[2]+v2[2];
    }
    template<class vettore3D>
    inline void diff(vettore3D& out, const vettore3D& v1,const vettore3D& v2){
        out[0] = v1[0]-v2[0];
        out[1] = v1[1]-v2[1];
        out[2] = v1[2]-v2[2];
    }
    template<class vettore3D>
    inline double dot(const vettore3D& v1,const vettore3D& v2){
        return (v1[0]*v2[0] +  v1[1]*v2[1]+ v1[2]*v2[2]);
    }
    template<class vettore3D>
    inline double norm(const vettore3D& v1){
            return std::sqrt(v1[0]*v1[0] +  v1[1]*v1[1]+ v1[2]*v1[2]);
    }
    template<class vettore3D>
    inline void cross(vettore3D& out, const vettore3D& v1,const vettore3D& v2){
        out[0] = v1[1]*v2[2] - v1[2]*v2[1];
        out[1] = v1[2]*v2[0] - v1[0]*v2[2];
        out[2] = v1[0]*v2[1] - v1[1]*v2[0];
    }
    template<class vettore3D>
    inline void vers(vettore3D& out, const vettore3D& in){
        double c = norm(in);
        for (int i = 0;i<3;++i) out[i] = in[i]/c;
    }
}

#endif // ARRAY3D_OPERATIONS_H
