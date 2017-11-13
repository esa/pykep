/*****************************************************************************
 *   Copyright (C) 2004-2018 The pykep development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *                                                                           *
 *   https://gitter.im/esa/pykep                                             *
 *   https://github.com/esa/pykep                                            *
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

#ifndef KEP_TOOLBOX_FB_PROP_H
#define KEP_TOOLBOX_FB_PROP_H

#include <cmath>

#include "../astro_constants.h"
#include "../core_functions/array3D_operations.h"

/// Propagate a fly-by hyperbola
/**
 * This function can be used to propagate a given fly-by hyperbola from its input conditions (entrance in the sphere of influence)
 * to its out conditions (exit from the sphere of influence)
 *
 * \param[out] v_out inertial velocity after encounter (cartesian)
 * \param[in]  v_in inertial velocity before encounter (cartesian)
 * \param[in]  v_pla inertial velocity of the planet at encounter (cartesian)
 * \param[in]  rp periplanet radius of the planetocentric hyperbola
 * \param[out] beta determines the hyperbola orbit plane orientation
 * \param[out] mu planet gravitationa parameter
 *
 * @see http://www.esa.int/gsp/ACT/doc/MAD/pub/ACT-RPR-MAD-2010-(CambridgePress)GlobalOptimizationAndSpacePruningForSpacecraftTrajectoryDesign.pdf
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

namespace kep_toolbox {

inline void fb_prop(array3D &v_out, const array3D& v_in, const array3D& v_pla, const double &rp, const double& beta, const double& mu)
{
    kep_toolbox::array3D v_rel_in = {{v_in[0]-v_pla[0], v_in[1]-v_pla[1], v_in[2]-v_pla[2]}};
    double v_rel_in2 = v_rel_in[0]*v_rel_in[0]+v_rel_in[1]*v_rel_in[1]+v_rel_in[2]*v_rel_in[2];
    double v_rel_in_norm = sqrt(v_rel_in2);
    double ecc = 1 + rp / mu * v_rel_in2;
    double delta = 2 * asin(1.0/ecc);
    kep_toolbox::array3D i_hat = {{v_rel_in[0] / v_rel_in_norm, v_rel_in[1] / v_rel_in_norm, v_rel_in[2] / v_rel_in_norm}};
    kep_toolbox::array3D j_hat = {{0,0,0}};
    kep_toolbox::cross(j_hat,i_hat,v_pla);
    kep_toolbox::vers(j_hat,j_hat);
    kep_toolbox::array3D k_hat= {{0,0,0}};
    kep_toolbox::cross(k_hat,i_hat,j_hat);
    v_out[0] = v_pla[0] + v_rel_in_norm * cos(delta) * i_hat[0] + v_rel_in_norm * cos(beta) * sin(delta) * j_hat[0] + v_rel_in_norm * sin(beta) * sin(delta) * k_hat[0];
    v_out[1] = v_pla[1] + v_rel_in_norm * cos(delta) * i_hat[1] + v_rel_in_norm * cos(beta) * sin(delta) * j_hat[1] + v_rel_in_norm * sin(beta) * sin(delta) * k_hat[1];
    v_out[2] = v_pla[2] + v_rel_in_norm * cos(delta) * i_hat[2] + v_rel_in_norm * cos(beta) * sin(delta) * j_hat[2] + v_rel_in_norm * sin(beta) * sin(delta) * k_hat[2];
    return;
}
} // namespace end

#endif // KEP_TOOLBOX_FB_PROP_H
