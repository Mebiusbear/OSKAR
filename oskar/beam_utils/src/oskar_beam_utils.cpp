/*
 * Copyright (c) 2022, The OSKAR Developers.
 * See the LICENSE file at the top-level directory of this distribution.
 */

#include <cmath>

#include "oskar_global.h"
#include "beam_utils/oskar_beam_utils.h"
#include "math/define_multiply.h"
#include "math/define_legendre_polynomial.h"
#include "telescope/station/element/define_evaluate_dipole_pattern.h"
#include "telescope/station/element/define_evaluate_spherical_wave.h"
#include "utility/oskar_kernel_macros.h"
#include "utility/oskar_vector_types.h"


OSKAR_EVALUATE_DIPOLE_PATTERN(evaluate_dipole_pattern, double, double2)
OSKAR_EVALUATE_SPHERICAL_WAVE_SUM(evaluate_spherical_wave_sum, double, double2, double4c)

namespace oskar{
namespace beam_utils{
void oskar_evaluate_dipole_pattern_double(
        const int num_points,
        const double* theta,
        const double* phi,
        const double freq_hz,
        const double dipole_length_m,
        std::complex<double>* pattern)
{
    const double kL = dipole_length_m * (M_PI * freq_hz / 299792458);
    const double cos_kL = cos(kL);
    const int stride = 4;
    double2* pattern_ptr = (double2*)pattern;
    evaluate_dipole_pattern(num_points, theta, phi, kL, cos_kL,
            stride, 0, 1, pattern_ptr, pattern_ptr);
}

void oskar_evaluate_spherical_wave_sum_double(
        double theta,
        double phi_x,
        double phi_y,
        int l_max,
        const std::complex<double>* alpha,
        std::complex<double>* pattern)
{
    const double4c* alpha_ptr = (double4c*)alpha;
    double4c* pattern_ptr = (double4c*)pattern;
    evaluate_spherical_wave_sum(1,
            &theta, &phi_x, &phi_y, l_max, alpha_ptr, 0, pattern_ptr);
}
} // namespace beam_utils
} // namespace oskar