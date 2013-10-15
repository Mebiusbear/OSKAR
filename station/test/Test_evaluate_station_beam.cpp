/*
 * Copyright (c) 2011-2013, The University of Oxford
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the University of Oxford nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <gtest/gtest.h>
#include <cuda_runtime_api.h>

#include <oskar_station.h>
#include <oskar_evaluate_station_beam_aperture_array.h>
#include <oskar_evaluate_station_beam_gaussian.h>
#include <oskar_evaluate_beam_horizontal_lmn.h>
#include <oskar_get_error_string.h>
#include <oskar_linspace.h>
#include <oskar_meshgrid.h>
#include <oskar_random_state.h>
#include <oskar_mem_binary_file_write.h>
#include <oskar_cuda_check_error.h>

#define TIMER_ENABLE 1
#include "utility/timer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

TEST(evaluate_station_beam, test_array_pattern)
{
    int error = 0;
    double gast = 0.0;

    // Construct a station model.
    double frequency = 30e6;
    int station_dim = 20;
    double station_size_m = 180.0;
    int num_antennas = station_dim * station_dim;
    oskar_Station* station = oskar_station_create(OSKAR_SINGLE,
            OSKAR_LOCATION_CPU, num_antennas, &error);
    oskar_station_resize_element_types(station, 1, &error);
    ASSERT_EQ(0, error) << oskar_get_error_string(error);

    // Set the station coordinates.
    oskar_station_set_position(station, 0.0, M_PI / 2.0, 0.0);
    float* x_pos = (float*) malloc(station_dim * sizeof(float));
    oskar_linspace_f(x_pos, -station_size_m/2.0, station_size_m/2.0, station_dim);
    oskar_meshgrid_f(
            oskar_mem_float(oskar_station_element_x_weights(station), &error),
            oskar_mem_float(oskar_station_element_y_weights(station), &error),
            x_pos, station_dim, x_pos, station_dim);
    free(x_pos);

    // Set the station beam direction.
    oskar_station_set_phase_centre(station,
            OSKAR_SPHERICAL_TYPE_EQUATORIAL, 0.0, M_PI / 2.0);

    // Set the station meta-data.
    oskar_Element* element = oskar_station_element(station, 0);
    oskar_element_set_element_type(element, OSKAR_ELEMENT_TYPE_ISOTROPIC);
    oskar_station_set_use_polarised_elements(station, 0);

    //error = oskar_station_save_configuration("temp_test_station.txt", &station_cpu);

    // Copy the station structure to the GPU and free the original structure.
    oskar_Station* station_gpu = oskar_station_create_copy(station,
            OSKAR_LOCATION_GPU, &error);
    oskar_station_free(station, &error);

    // Evaluate horizontal l,m positions at which to generate the beam pattern.
    int image_size = 301;
    double fov_deg = 30.0;
    int num_pixels = image_size * image_size;

    // Generate horizontal lm coordinates for the beam pattern.
    oskar_Mem beam_pattern, h_l, h_m, h_n, d_l, d_m, d_n;
    oskar_mem_init(&h_l, OSKAR_SINGLE, OSKAR_LOCATION_CPU, num_pixels, 1, &error);
    oskar_mem_init(&h_m, OSKAR_SINGLE, OSKAR_LOCATION_CPU, num_pixels, 1, &error);
    oskar_mem_init(&h_n, OSKAR_SINGLE, OSKAR_LOCATION_CPU, num_pixels, 1, &error);
    float* lm = (float*)malloc(image_size * sizeof(float));
    double lm_max = sin(fov_deg * M_PI / 180.0);
    oskar_linspace_f(lm, -lm_max, lm_max, image_size);
    oskar_meshgrid_f(oskar_mem_float(&h_l, &error),
            oskar_mem_float(&h_m, &error), lm, image_size, lm, image_size);
    free(lm);

    // Copy horizontal lm coordinates to GPU.
    oskar_mem_init_copy(&d_l, &h_l, OSKAR_LOCATION_GPU, &error);
    oskar_mem_init_copy(&d_m, &h_m, OSKAR_LOCATION_GPU, &error);
    oskar_mem_init_copy(&d_n, &h_n, OSKAR_LOCATION_GPU, &error);

    // Initialise the random number generator.
    int seed = 0;
    oskar_RandomState* random_state = oskar_random_state_create(num_antennas,
            seed, 0, 0, &error);
    ASSERT_EQ(0, error) << oskar_get_error_string(error);

    // Allocate weights work array.
    oskar_StationWork* work = oskar_station_work_create(OSKAR_SINGLE,
            OSKAR_LOCATION_GPU, &error);

    // Create memory for the beam pattern.
    oskar_mem_init(&beam_pattern, OSKAR_SINGLE_COMPLEX, OSKAR_LOCATION_GPU,
            num_pixels, 1, &error);

    ASSERT_EQ(0, oskar_station_array_is_3d(station_gpu));
    TIMER_START
    oskar_evaluate_station_beam_aperture_array(&beam_pattern, station_gpu,
            num_pixels, &d_l, &d_m, &d_n, gast, frequency, work,
            random_state, &error);
    cudaDeviceSynchronize();
    TIMER_STOP("Finished aperture array station beam (2D)");
    ASSERT_EQ(0, error) << oskar_get_error_string(error);
    oskar_station_set_element_coords(station_gpu, 0, 0., 0., 1., 0., 0., 0., &error);
    oskar_station_set_element_coords(station_gpu, 0, 0., 0., 0., 0., 0., 0., &error);
    ASSERT_EQ(1, oskar_station_array_is_3d(station_gpu));
    TIMER_START
    oskar_evaluate_station_beam_aperture_array(&beam_pattern, station_gpu,
            num_pixels, &d_l, &d_m, &d_n, gast, frequency, work,
            random_state, &error);
    cudaDeviceSynchronize();
    TIMER_STOP("Finished aperture array station beam (3D)");
    ASSERT_EQ(0, error) << oskar_get_error_string(error);
    oskar_station_free(station_gpu, &error);

    // Save beam to file for plotting.
    const char* filename = "temp_test_beam_pattern.txt";
    FILE* file = fopen(filename, "w");
    oskar_mem_write_ascii(file, 3, num_pixels, &error,
            &h_l, &h_m, &beam_pattern);
    fclose(file);
    remove(filename);

    /*--------------------------------------------------------------------------
        data = dlmread('temp_test_beam_pattern.txt');
        imagesc(log10(reshape(data(:,3), 301, 301).^2));
    --------------------------------------------------------------------------*/
    oskar_random_state_free(random_state, &error);
    oskar_station_work_free(work, &error);
    oskar_mem_free(&beam_pattern, &error);
    oskar_mem_free(&h_l, &error);
    oskar_mem_free(&h_m, &error);
    oskar_mem_free(&h_n, &error);
    oskar_mem_free(&d_l, &error);
    oskar_mem_free(&d_m, &error);
    oskar_mem_free(&d_n, &error);
    ASSERT_EQ(0, error) << oskar_get_error_string(error);
}


TEST(evaluate_station_beam, gaussian)
{
    int error = 0;
    double fwhm = 1.0;
    int size = 256;
    int num_points = size * size;
    double lm_minmax = 0.2;
    bool save_results = false;

    // Double CPU
    {
        int type = OSKAR_DOUBLE;
        int location = OSKAR_LOCATION_CPU;
        double* x = (double*)malloc(size * sizeof(double));
        oskar_linspace_d(x, -lm_minmax, lm_minmax, size);
        oskar_Mem l, m, beam, horizon_mask;
        oskar_mem_init(&l, type, location, num_points, 1, &error);
        oskar_mem_init(&m, type, location, num_points, 1, &error);
        oskar_mem_init(&horizon_mask, type, location, num_points, 1, &error);
        oskar_mem_init(&beam, type | OSKAR_COMPLEX, location, num_points, 1, &error);
        oskar_meshgrid_d(oskar_mem_double(&l, &error),
                oskar_mem_double(&m, &error), x, size, x, size);
        free(x);

        oskar_evaluate_station_beam_gaussian(&beam, num_points, &l, &m,
                &horizon_mask, fwhm, &error);
        ASSERT_EQ(0, error) << oskar_get_error_string(error);

        // Write output to file.
        if (save_results)
        {
            const char* filename = "temp_test_beam_double_cpu.dat";
            remove(filename);
            oskar_mem_binary_file_write(&beam, filename, 0, 0, 0,
                    (int)oskar_mem_length(&beam), &error);
            ASSERT_EQ(0, error) << oskar_get_error_string(error);
        }
        oskar_mem_free(&l, &error);
        oskar_mem_free(&m, &error);
        oskar_mem_free(&horizon_mask, &error);
        oskar_mem_free(&beam, &error);
    }

    // Single CPU
    {
        int type = OSKAR_SINGLE;
        int location = OSKAR_LOCATION_CPU;
        float* x = (float*)malloc(size * sizeof(float));
        oskar_linspace_f(x, -lm_minmax, lm_minmax, size);
        oskar_Mem l, m, beam, horizon_mask;
        oskar_mem_init(&l, type, location, num_points, 1, &error);
        oskar_mem_init(&m, type, location, num_points, 1, &error);
        oskar_mem_init(&horizon_mask, type, location, num_points, 1, &error);
        oskar_mem_init(&beam, type | OSKAR_COMPLEX, location, num_points, 1, &error);
        oskar_meshgrid_f(oskar_mem_float(&l, &error),
                oskar_mem_float(&m, &error), x, size, x, size);
        free(x);

        oskar_evaluate_station_beam_gaussian(&beam, num_points, &l, &m,
                &horizon_mask, fwhm, &error);
        ASSERT_EQ(0, error) << oskar_get_error_string(error);

        // Write output to file.
        if (save_results)
        {
            const char* filename = "temp_test_beam_single_cpu.dat";
            remove(filename);
            oskar_mem_binary_file_write(&beam, filename, 0, 0, 0,
                    (int)oskar_mem_length(&beam), &error);
            ASSERT_EQ(0, error) << oskar_get_error_string(error);
        }
        oskar_mem_free(&l, &error);
        oskar_mem_free(&m, &error);
        oskar_mem_free(&horizon_mask, &error);
        oskar_mem_free(&beam, &error);
    }

    // Double GPU
    {
        int type = OSKAR_DOUBLE;
        int location = OSKAR_LOCATION_GPU;
        double* x = (double*)malloc(size * sizeof(double));
        oskar_linspace_d(x, -lm_minmax, lm_minmax, size);
        oskar_Mem h_l, h_m, l, m, horizon_mask, beam;
        oskar_mem_init(&h_l, type, OSKAR_LOCATION_CPU, num_points, 1, &error);
        oskar_mem_init(&h_m, type, OSKAR_LOCATION_CPU, num_points, 1, &error);
        oskar_mem_init(&horizon_mask, type, location, num_points, 1, &error);
        oskar_mem_init(&beam, type | OSKAR_COMPLEX, location, num_points, 1, &error);
        oskar_meshgrid_d(oskar_mem_double(&h_l, &error),
                oskar_mem_double(&h_m, &error), x, size, x, size);
        free(x);

        oskar_mem_init_copy(&l, &h_l, location, &error);
        oskar_mem_init_copy(&m, &h_m, location, &error);
        oskar_evaluate_station_beam_gaussian(&beam, num_points, &l, &m,
                &horizon_mask, fwhm, &error);
        ASSERT_EQ(0, error) << oskar_get_error_string(error);

        // Write output to file.
        if (save_results)
        {
            const char* filename = "temp_test_beam_double_gpu.dat";
            remove(filename);
            oskar_mem_binary_file_write(&beam, filename, 0, 0, 0,
                    (int)oskar_mem_length(&beam), &error);
            ASSERT_EQ(0, error) << oskar_get_error_string(error);
        }
        oskar_mem_free(&h_l, &error);
        oskar_mem_free(&h_m, &error);
        oskar_mem_free(&l, &error);
        oskar_mem_free(&m, &error);
        oskar_mem_free(&horizon_mask, &error);
        oskar_mem_free(&beam, &error);
    }

    // Single GPU
    {
        int type = OSKAR_SINGLE;
        int location = OSKAR_LOCATION_GPU;
        float* x = (float*)malloc(size * sizeof(float));
        oskar_linspace_f(x, -lm_minmax, lm_minmax, size);
        oskar_Mem h_l, h_m, l, m, horizon_mask, beam;
        oskar_mem_init(&h_l, type, OSKAR_LOCATION_CPU, num_points, 1, &error);
        oskar_mem_init(&h_m, type, OSKAR_LOCATION_CPU, num_points, 1, &error);
        oskar_mem_init(&horizon_mask, type, location, num_points, 1, &error);
        oskar_mem_init(&beam, type | OSKAR_COMPLEX, location, num_points, 1, &error);
        oskar_meshgrid_f(oskar_mem_float(&h_l, &error),
                oskar_mem_float(&h_m, &error), x, size, x, size);
        free(x);

        oskar_mem_init_copy(&l, &h_l, location, &error);
        oskar_mem_init_copy(&m, &h_m, location, &error);
        oskar_evaluate_station_beam_gaussian(&beam, num_points, &l, &m,
                &horizon_mask, fwhm, &error);
        ASSERT_EQ(0, error) << oskar_get_error_string(error);

        // Write output to file.
        if (save_results)
        {
            const char* filename = "temp_test_beam_single_gpu.dat";
            remove(filename);
            oskar_mem_binary_file_write(&beam, filename, 0, 0, 0,
                    (int)oskar_mem_length(&beam), &error);
            ASSERT_EQ(0, error) << oskar_get_error_string(error);
        }
        oskar_mem_free(&h_l, &error);
        oskar_mem_free(&h_m, &error);
        oskar_mem_free(&l, &error);
        oskar_mem_free(&m, &error);
        oskar_mem_free(&horizon_mask, &error);
        oskar_mem_free(&beam, &error);
    }
}